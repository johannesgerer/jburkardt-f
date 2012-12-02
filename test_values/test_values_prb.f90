program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_VALUES_PRB.
!
!  Discussion:
!
!    TEST_VALUES_PRB calls the TEST_VALUE routines.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_VALUES_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_VALUES library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test0035 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test008 ( )
  call test009 ( )
  call test0093 ( )
  call test0095 ( )

  call test010 ( )
  call test011 ( )
  call test0114 ( )
  call test01145 ( )
  call test0115 ( )
  call test01155 ( )
  call test0116 ( )
  call test012 ( )
  call test0123 ( )
  call test0127 ( )
  call test0128 ( )
  call test013 ( )
  call test0134 ( )
  call test0135 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )
  call test017 ( )
  call test018 ( )
  call test0185 ( )
  call test019 ( )
  call test0195 ( )

  call test020 ( )
  call test0205 ( )
  call test021 ( )
  call test022 ( )
  call test023 ( )
  call test024 ( )
  call test025 ( )
  call test026 ( )
  call test0265 ( )
  call test027 ( )
  call test028 ( )
  call test029 ( )

  call test030 ( )
  call test0305 ( )
  call test031 ( )
  call test032 ( )
  call test033 ( )
  call test034 ( )
  call test035 ( )
  call test036 ( )
  call test0365 ( )
  call test037 ( )
  call test038 ( )
  call test039 ( )
  call test0395 ( )

  call test040 ( )
  call test041 ( )
  call test042 ( )
  call test0425 ( )
  call test043 ( )
  call test044 ( )
  call test0445 ( )
  call test045 ( )
  call test046 ( )
  call test0465 ( )
  call test047 ( )
  call test048 ( )
  call test049 ( )

  call test050 ( )
  call test051 ( )
  call test05125 ( )
  call test0515 ( )
  call test0517 ( )
  call test0519 ( )
  call test052 ( )
  call test053 ( )
  call test054 ( )
  call test055 ( )
  call test056 ( )
  call test057 ( )
  call test0575 ( )
  call test058 ( )
  call test059 ( )

  call test060 ( )
  call test061 ( )
  call test062 ( )
  call test063 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test0665 ( )
  call test067 ( )
  call test068 ( )
  call test0685 ( )
  call test069 ( )

  call test070 ( )
  call test071 ( )
  call test072 ( )
  call test073 ( )
  call test074 ( )
  call test075 ( )
  call test0755 ( )
  call test0756 ( )
  call test076 ( )
  call test077 ( )
  call test078 ( )
  call test079 ( )

  call test080 ( )
  call test081 ( )
  call test082 ( )
  call test083 ( )
  call test0835 ( )
  call test084 ( )
  call test0843 ( )
  call test0845 ( )
  call test085 ( )
  call test0855 ( )
  call test086 ( )
  call test087 ( )
  call test088 ( )
  call test089 ( )

  call test090 ( )
  call test091 ( )
  call test092 ( )
  call test093 ( )
  call test094 ( )
  call test0945 ( )
  call test095 ( )
  call test096 ( )
  call test097 ( )
  call test0972 ( )
  call test0973 ( )
  call test0974 ( )
  call test0975 ( )
  call test098 ( )
  call test099 ( )
  call test0995 ( )

  call test100 ( )
  call test101 ( )
  call test1015 ( )
  call test1016 ( )
  call test102 ( )
  call test103 ( )
  call test1035 ( )
  call test104 ( )
  call test1037 ( )
  call test105 ( )
  call test106 ( )
  call test107 ( )
  call test108 ( )
  call test10825 ( )
  call test10850 ( )
  call test10875 ( )
  call test109 ( )

  call test110 ( )
  call test1105 ( )
  call test111 ( )
  call test112 ( )
  call test113 ( )
  call test1135 ( )
  call test114 ( )
  call test115 ( )
  call test116 ( )
  call test117 ( )
  call test118 ( )
  call test1185 ( )
  call test119 ( )

  call test120 ( )
  call test121 ( )
  call test122 ( )
  call test123 ( )
  call test124 ( )
  call test125 ( )
  call test1255 ( )
  call test126 ( )
  call test127 ( )
  call test1275 ( )
  call test128 ( )
  call test1283 ( )
  call test1285 ( )
  call test129 ( )

  call test131 ( )
  call test132 ( )
  call test1325 ( )
  call test130 ( )
  call test133 ( )
  call test134 ( )
  call test135 ( )
  call test136 ( )
  call test137 ( )
  call test138 ( )
  call test139 ( )

  call test140 ( )
  call test141 ( )
  call test1415 ( )
  call test142 ( )
  call test143 ( )
  call test144 ( )
  call test1445 ( )
  call test1447 ( )
  call test145 ( )
  call test146 ( )
  call test1465 ( )
  call test147 ( )
  call test148 ( )
  call test149 ( )

  call test150 ( )
  call test151 ( )
  call test152 ( )
  call test153 ( )
  call test154 ( )
  call test1545 ( )
  call test155 ( )
  call test156 ( )
  call test157 ( )
  call test1575 ( )
  call test158 ( )
  call test159 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_VALUES_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests ABRAM0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001:'
  write ( *, '(a)' ) '  ABRAM0_VALUES returns values of '
  write ( *, '(a)' ) '  the Abramowitz function of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Abram0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests ABRAM1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002:'
  write ( *, '(a)' ) '  ABRAM1_VALUES returns values of '
  write ( *, '(a)' ) '  the Abramowitz function of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Abram1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests ABRAM2_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003:'
  write ( *, '(a)' ) '  ABRAM2_VALUES returns values of '
  write ( *, '(a)' ) '  the Abramowitz function of order 2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Abram2'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram2_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0035 ( )

!*****************************************************************************80
!
!! TEST0035 tests AGM_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0035:'
  write ( *, '(a)' ) '  AGM_VALUES returns values of '
  write ( *, '(a)' ) '  the arithmetic geometric mean function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          A               B            AGM(A,B)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call agm_values ( n_data, a, b, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) a, b, fx

  end do

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests AIRY_AI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ai
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004:'
  write ( *, '(a)' ) '  AIRY_AI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Ai(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Ai'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_values ( n_data, x, ai )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, ai

  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests AIRY_AI_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005:'
  write ( *, '(a)' ) '  AIRY_AI_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Airy function Ai(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Ai_Int'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests AIRY_AI_PRIME_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) aip
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006:'
  write ( *, '(a)' ) '  AIRY_AI_PRIME_VALUES returns values of '
  write ( *, '(a)' ) '  the derivative of the Airy functions A''(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           AiP'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_prime_values ( n_data, x, aip )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, aip

  end do

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests AIRY_BI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bi
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007:'
  write ( *, '(a)' ) '  AIRY_BI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Bi(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Bi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_values ( n_data, x, bi )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, bi

  end do

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests AIRY_BI_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008:'
  write ( *, '(a)' ) '  AIRY_BI_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Airy function Bi(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Bi_Int'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests AIRY_BI_PRIME_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bip
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009:'
  write ( *, '(a)' ) '  AIRY_BI_PRIME_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function derivative B''(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BiP'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_prime_values ( n_data, x, bip )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, bip

  end do

  return
end
subroutine test0093 ( )

!*****************************************************************************80
!
!! TEST0093 tests AIRY_CAI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 8 ) ai
  integer ( kind = 4 ) n_data
  complex ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004:'
  write ( *, '(a)' ) '  AIRY_CAI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Ai(X) with complex argument'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X                         Ai'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_cai_values ( n_data, x, ai )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,(f14.6,f14.6),2x,(g24.16,g24.16))' ) x, ai

  end do

  return
end
subroutine test0095 ( )

!*****************************************************************************80
!
!! TEST0095 tests AIRY_CBI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  complex ( kind = 8 ) bi
  integer ( kind = 4 ) n_data
  complex ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005:'
  write ( *, '(a)' ) '  AIRY_CBI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Bi(X) with complex argument'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X                         Bi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_cbi_values ( n_data, x, bi )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,(f14.6,f14.6),2x,(g24.16,g24.16))' ) x, bi

  end do

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests AIRY_GI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010:'
  write ( *, '(a)' ) '  AIRY_GI_VALUES returns values of '
  write ( *, '(a)' ) '  the modified Airy function Gi(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Gi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_gi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests AIRY_HI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011:'
  write ( *, '(a)' ) '  AIRY_HI_VALUES returns values of '
  write ( *, '(a)' ) '  the modified Airy function Hi(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Hi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_hi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0114 ( )

!*****************************************************************************80
!
!! TEST0114 tests ARCCOS_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0114:'
  write ( *, '(a)' ) '  ARCCOS_VALUES returns values of '
  write ( *, '(a)' ) '  the arc cosine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arccos_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test01145 ( )

!*****************************************************************************80
!
!! TEST01145 tests ARCCOSH_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01145:'
  write ( *, '(a)' ) '  ARCCOSH_VALUES returns values of '
  write ( *, '(a)' ) '  the hyperbolic arc cosine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arccosh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0115 ( )

!*****************************************************************************80
!
!! TEST0115 tests ARCSIN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0115:'
  write ( *, '(a)' ) '  ARCSIN_VALUES returns values of '
  write ( *, '(a)' ) '  the arc sine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arcsin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test01155 ( )

!*****************************************************************************80
!
!! TEST01155 tests ARCSINH_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01155:'
  write ( *, '(a)' ) '  ARCSINH_VALUES returns values of '
  write ( *, '(a)' ) '  the hyperbolic arc sine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arcsinh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0116 ( )

!*****************************************************************************80
!
!! TEST0116 tests ARCTAN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0116:'
  write ( *, '(a)' ) '  ARCTAN_VALUES returns values of '
  write ( *, '(a)' ) '  the arc tangent'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctan_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests ARCTAN_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012:'
  write ( *, '(a)' ) '  ARCTAN_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the arctangent integral'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctan_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0123 ( )

!*****************************************************************************80
!
!! TEST0123 tests ARCTANH_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0123:'
  write ( *, '(a)' ) '  ARCTANH_VALUES returns values of '
  write ( *, '(a)' ) '  the hyperbolic arc tangent'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctanh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0127 ( )

!*****************************************************************************80
!
!! TEST0127 tests BEI0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0127:'
  write ( *, '(a)' ) '  BEI0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BEI of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BEI0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bei0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0128 ( )

!*****************************************************************************80
!
!! TEST0128 tests BEI1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0128:'
  write ( *, '(a)' ) '  BEI1_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BEI of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BEI1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bei1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests BELL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013:'
  write ( *, '(a)' ) '  BELL_VALUES returns values of '
  write ( *, '(a)' ) '  the Bell numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        BELL(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bell_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, c

  end do

  return
end
subroutine test0134 ( )

!*****************************************************************************80
!
!! TEST0134 tests BER0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0134:'
  write ( *, '(a)' ) '  BER0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BER of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BEI0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ber0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0135 ( )

!*****************************************************************************80
!
!! TEST0135 tests BER1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0135:'
  write ( *, '(a)' ) '  BER1_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BER of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BER1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ber1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests BERNOULLI_NUMBER_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014:'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER_VALUES returns values of '
  write ( *, '(a)' ) '  the Bernoulli numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        BERNOULLI_NUMBER(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,6x,g24.16)' ) n, c

  end do

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests BERNOULLI_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015:'
  write ( *, '(a)' ) '  BERNOULLI_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Bernoulli polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            BERNOULLI_POLY(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernoulli_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests BERNSTEIN_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016:'
  write ( *, '(a)' ) '  BERNSTEIN_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Bernstein polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         K          X           BERNSTEIN(N,K)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernstein_poly_values ( n_data, n, k, x, b )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, k, x, b

  end do

  return
end
subroutine test017 ( )

!*****************************************************************************80
!
!! TEST017 tests BESSEL_I0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST017:'
  write ( *, '(a)' ) '  BESSEL_I0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel I0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            I0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test018 ( )

!*****************************************************************************80
!
!! TEST018 tests BESSEL_I0_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST018:'
  write ( *, '(a)' ) '  BESSEL_I0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel I0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            I0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0185 ( )

!*****************************************************************************80
!
!! TEST0185 tests BESSEL_I0_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0185:'
  write ( *, '(a)' ) '  BESSEL_I0_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel i0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           i0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test019 ( )

!*****************************************************************************80
!
!! TEST019 tests BESSEL_I1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST019:'
  write ( *, '(a)' ) '  BESSEL_I1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel I1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            I1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0195 ( )

!*****************************************************************************80
!
!! TEST0195 tests BESSEL_I1_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0195:'
  write ( *, '(a)' ) '  BESSEL_I1_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel i1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           i1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i1_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests BESSEL_IN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020:'
  write ( *, '(a)' ) '  BESSEL_IN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel In function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X         I(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_in_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0205 ( )

!*****************************************************************************80
!
!! TEST0205 tests BESSEL_IX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0205:'
  write ( *, '(a)' ) '  BESSEL_IX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel In function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         I(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_ix_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests BESSEL_J0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021:'
  write ( *, '(a)' ) '  BESSEL_J0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel J0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           J0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests BESSEL_J0_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022:'
  write ( *, '(a)' ) '  BESSEL_J0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel J0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           J0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests BESSEL_J0_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023:'
  write ( *, '(a)' ) '  BESSEL_J0_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel j0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           j0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests BESSEL_J1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024:'
  write ( *, '(a)' ) '  BESSEL_J1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel J1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           J1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests BESSEL_J1_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025:'
  write ( *, '(a)' ) '  BESSEL_J1_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel j1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           j1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j1_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests BESSEL_JN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026:'
  write ( *, '(a)' ) '  BESSEL_JN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Jn function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X           J(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0265 ( )

!*****************************************************************************80
!
!! TEST0265 tests BESSEL_JX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0265:'
  write ( *, '(a)' ) '  BESSEL_JX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Jn function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         J(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jx_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests BESSEL_K0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027:'
  write ( *, '(a)' ) '  BESSEL_K0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel K0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            K0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests BESSEL_K0_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028:'
  write ( *, '(a)' ) '  BESSEL_K0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel K0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            K0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests BESSEL_K1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029:'
  write ( *, '(a)' ) '  BESSEL_K1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel K1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            K1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test030 ( )

!*****************************************************************************80
!
!! TEST030 tests BESSEL_KN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TES030:'
  write ( *, '(a)' ) '  BESSEL_KN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Kn function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            K(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_kn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0305 ( )

!*****************************************************************************80
!
!! TEST0305 tests BESSEL_KX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0305:'
  write ( *, '(a)' ) '  BESSEL_KX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Kn function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         K(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_kx_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests BESSEL_Y0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031:'
  write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Y0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            X            Y0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests BESSEL_Y0_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032:'
  write ( *, '(a)' ) '  BESSEL_Y0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel Y0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Y0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test033 ( )

!*****************************************************************************80
!
!! TEST033 tests BESSEL_Y0_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033:'
  write ( *, '(a)' ) '  BESSEL_Y0_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel Y0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            y0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests BESSEL_Y1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034:'
  write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Y1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Y1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! TEST035 tests BESSEL_Y1_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035:'
  write ( *, '(a)' ) '  BESSEL_Y1_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel Y1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            y1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y1_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests BESSEL_YN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036:'
  write ( *, '(a)' ) '  BESSEL_YN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Yn function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X              Y(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_yn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0365 ( )

!*****************************************************************************80
!
!! TEST0365 tests BESSEL_YX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0365:'
  write ( *, '(a)' ) '  BESSEL_YX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Yn function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         Y(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_yx_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 tests BETA_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037:'
  write ( *, '(a)' ) '  BETA_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Beta CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B               X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g14.6)' ) a, b, x, fx

  end do

  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests BETA_INC_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038:'
  write ( *, '(a)' ) '  BETA_INC_VALUES returns values of '
  write ( *, '(a)' ) '  the incomplete Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B               X           BETA_INC(A,B)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g14.6)' ) a, b, x, fx

  end do

  return
end
subroutine test039 ( )

!*****************************************************************************80
!
!! TEST039 tests BETA_LOG_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fxy
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039:'
  write ( *, '(a)' ) '  BETA_LOG_VALUES returns values of '
  write ( *, '(a)' ) '  the logarithm of the Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y           BETA_LOG(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_log_values ( n_data, x, y, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) x, y, fxy

  end do

  return
end
subroutine test0395 ( )

!*****************************************************************************80
!
!! TEST0395 tests BETA_NONCENTRAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0395:'
  write ( *, '(a)' ) '  BETA_NONCENTRAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Beta CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '        A         B     LAMBDA          X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.2,2x,f8.2,2x,f8.2,2x,f14.6,2x,g24.16)' ) &
    a, b, lambda, x, fx

  end do

  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests BETA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fxy
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040:'
  write ( *, '(a)' ) '  BETA_VALUES returns values of '
  write ( *, '(a)' ) '  the Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y            BETA(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_values ( n_data, x, y, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) x, y, fxy

  end do

  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 tests BINOMIAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041:'
  write ( *, '(a)' ) '  BINOMIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the binomial numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A         B        C(A,B)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_values ( n_data, a, b, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) a, b, c

  end do

  return
end
subroutine test042 ( )

!*****************************************************************************80
!
!! TEST042 tests BINOMIAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042:'
  write ( *, '(a)' ) '  BINOMIAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Binomial Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A        B            X       CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f10.4,2x,i8,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine test0425 ( )

!*****************************************************************************80
!
!! TEST0425 tests BIVARIATE_NORMAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fxy
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0425:'
  write ( *, '(a)' ) '  BIVARIATE_NORMAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the bivariate normal CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X               Y               R           F(R)(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bivariate_normal_cdf_values ( n_data, x, y, r, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g14.6)' ) x, y, r, fxy

  end do

  return
end
subroutine test043 ( )

!*****************************************************************************80
!
!! TEST043 tests CATALAN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043:'
  write ( *, '(a)' ) '  CATALAN_VALUES returns values of '
  write ( *, '(a)' ) '  the Catalan numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          C(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call catalan_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, c

  end do

  return
end
subroutine test044 ( )

!*****************************************************************************80
!
!! TEST044 tests CAUCHY_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST044:'
  write ( *, '(a)' ) '  CAUCHY_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Cauchy Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       MU              SIGMA           X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cauchy_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine test0445 ( )

!*****************************************************************************80
!
!! TEST0445 tests CBRT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0445:'
  write ( *, '(a)' ) '  CBRT_VALUES returns values of the cube root.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X         CBRT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cbrt_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests CHEBY_T_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045:'
  write ( *, '(a)' ) '  CHEBY_T_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Chebyshev T polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         T(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_t_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test046 ( )

!*****************************************************************************80
!
!! TEST046 tests CHEBY_U_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046:'
  write ( *, '(a)' ) '  CHEBY_U_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Chebyshev U polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X          U(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_u_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0465 ( )

!*****************************************************************************80
!
!! TEST0465 tests CHI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0465:'
  write ( *, '(a)' ) '  CHI_VALUES returns values of '
  write ( *, '(a)' ) '  the Hyperbolic Cosine Integral function CHI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           CHI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test047 ( )

!*****************************************************************************80
!
!! TEST047 tests CHI_SQUARE_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047:'
  write ( *, '(a)' ) '  CHI_SQUARE_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Chi-Squared Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine test048 ( )

!*****************************************************************************80
!
!! TEST048 tests CHI_SQUARE_NONCENTRAL_CDF_VALUES.
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

  integer ( kind = 4 ) df
  real ( kind = 8 ) fx
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048:'
  write ( *, '(a)' ) '  CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the noncentral Chi-Squared Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        DF        LAMBDA        X         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_noncentral_cdf_values ( n_data, df, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i10,2x,f10.4,2x,f10.4,2x,g24.16)' ) df, lambda, x, fx

  end do

  return
end
subroutine test049 ( )

!*****************************************************************************80
!
!! TEST049 tests CI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049:'
  write ( *, '(a)' ) '  CI_VALUES returns values of '
  write ( *, '(a)' ) '  the Cosine Integral function CI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           CI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ci_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test050 ( )

!*****************************************************************************80
!
!! TEST050 tests CIN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050:'
  write ( *, '(a)' ) '  CIN_VALUES returns values of '
  write ( *, '(a)' ) '  the Cosine Integral function CIN(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X            CIN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test051 ( )

!*****************************************************************************80
!
!! TEST051 tests CLAUSEN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051:'
  write ( *, '(a)' ) '  CLAUSEN_VALUES returns values of '
  write ( *, '(a)' ) '  Clausen''s integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            FX'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call clausen_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test05125 ( )

!*****************************************************************************80
!
!! TEST05125 tests CLEBSCH_GORDAN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) j1
  real ( kind = 8 ) j2
  real ( kind = 8 ) j3
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) m3
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05125:'
  write ( *, '(a)' ) '  CLEBSCH_GORDAN_VALUES returns values of '
  write ( *, '(a)' ) '  the Clebsch-Gordan coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      M1      M2      M3      CG'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call clebsch_gordan_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, m1, m2, m3, fx

  end do

  return
end
subroutine test0515 ( )

!*****************************************************************************80
!
!! TEST0515 tests COLLATZ_COUNT_VALUES.
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

  integer ( kind = 4 ) count
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0515:'
  write ( *, '(a)' ) '  COLLATZ_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the length of the Collatz sequence that'
  write ( *, '(a)' ) '  starts at N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N      COLLATZ_COUNT(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call collatz_count_values ( n_data, n, count )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, count

  end do

  return
end
subroutine test0517 ( )

!*****************************************************************************80
!
!! TEST0517 tests COS_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0517:'
  write ( *, '(a)' ) '  COS_VALUES returns values of the cosine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           COS(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cos_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0519 ( )

!*****************************************************************************80
!
!! TEST0519 tests COSH_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0519:'
  write ( *, '(a)' ) &
    '  COSH_VALUES returns values of the hyperbolic cosine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           COSH(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cosh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests CP_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cp
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052:'
  write ( *, '(a)' ) '  CP_VALUES returns values of '
  write ( *, '(a)' ) '  the specific heat CP '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          T               P          CP(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cp_values ( n_data, tc, p, cp )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, cp

  end do

  return
end
subroutine test053 ( )

!*****************************************************************************80
!
!! TEST053 tests DAWSON_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST053:'
  write ( *, '(a)' ) '  DAWSON_VALUES returns values of '
  write ( *, '(a)' ) '  Dawson''s Integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            DAWSON(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dawson_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests DEBYE1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054:'
  write ( *, '(a)' ) '  DEBYE1_VALUES returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests DEBYE2_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055:'
  write ( *, '(a)' ) '  DEBYE2_VALUES returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE2(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye2_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test056 ( )

!*****************************************************************************80
!
!! TEST056 tests DEBYE3_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST056:'
  write ( *, '(a)' ) '  DEBYE3_VALUES returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE3(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye3_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test057 ( )

!*****************************************************************************80
!
!! TEST057 tests DEBYE4_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057:'
  write ( *, '(a)' ) '  DEBYE4_VALUES returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE4(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye4_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0575 ( )

!*****************************************************************************80
!
!! TEST0575 tests DILOGARITHM_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) p
  integer ( kind = 4 ) q

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0575:'
  write ( *, '(a)' ) '  DEDEKIND_SUM_VALUES returns values of the'
  write ( *, '(a)' ) '  Dedekind sum function: (N/D) = Dedekind_Sum ( P, Q ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P     Q     N     D'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dedekind_sum_values ( n_data, p, q, n, d )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) p, q, n, d

  end do

  return
end
subroutine test058 ( )

!*****************************************************************************80
!
!! TEST058 tests DIELECTRIC_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) eps
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058:'
  write ( *, '(a)' ) '  DIELECTRIC_VALUES returns values of '
  write ( *, '(a)' ) '  the dielectric function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          T            P            EPS(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dielectric_values ( n_data, tc, p, eps )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, eps

  end do

  return
end
subroutine test059 ( )

!*****************************************************************************80
!
!! TEST059 tests DILOGARITHM_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059:'
  write ( *, '(a)' ) '  DILOGARITHM_VALUES returns values of'
  write ( *, '(a)' ) '  the dilogarithm function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          DILOGARITHM(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dilogarithm_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests E1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060:'
  write ( *, '(a)' ) '  E1_VALUES returns values of'
  write ( *, '(a)' ) '  the exponential integral function E1(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          E1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call e1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test061 ( )

!*****************************************************************************80
!
!! TEST061 tests EI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061:'
  write ( *, '(a)' ) '  EI_VALUES returns values of'
  write ( *, '(a)' ) '  the exponential integral function EI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          EI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ei_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests ELLIPTIC_EA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062:'
  write ( *, '(a)' ) '  ELLIPTIC_EA_VALUES returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameter angle ALPHA in degrees.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          ALPHA        EA(ALPHA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_ea_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test063 ( )

!*****************************************************************************80
!
!! TEST063 tests ELLIPTIC_EM_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063:'
  write ( *, '(a)' ) '  ELLIPTIC_EM_VALUES returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameter modulus M.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          M            EM(M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_em_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test064 ( )

!*****************************************************************************80
!
!! TEST064 tests ELLIPTIC_KA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064:'
  write ( *, '(a)' ) '  ELLIPTIC_KA_VALUES returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameter angle ALPHA in degrees.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          ALPHA        KA(ALPHA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_ka_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests ELLIPTIC_KM_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065:'
  write ( *, '(a)' ) '  ELLIPTIC_KM_VALUES returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameter modulus M.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          M            KM(M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_km_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test066 ( )

!*****************************************************************************80
!
!! TEST066 tests ERF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066:'
  write ( *, '(a)' ) '  ERF_VALUES returns values of'
  write ( *, '(a)' ) '  the error function ERF(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           ERF(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0665 ( )

!*****************************************************************************80
!
!! TEST0665 tests ERFC_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0665:'
  write ( *, '(a)' ) '  ERFC_VALUES returns values of'
  write ( *, '(a)' ) '  the complementary error function ERFC(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           ERFC(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erfc_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test067 ( )

!*****************************************************************************80
!
!! TEST067 tests EULER_NUMBER_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) c
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067:'
  write ( *, '(a)' ) '  EULER_NUMBER_VALUES returns values of '
  write ( *, '(a)' ) '  the Euler numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N        EULER_NUMBER(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_number_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, c

  end do

  return
end
subroutine test068 ( )

!*****************************************************************************80
!
!! TEST068 tests EULER_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068:'
  write ( *, '(a)' ) '  EULER_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Euler polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X             EULER_POLY(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0685 ( )

!*****************************************************************************80
!
!! TEST0685 tests EXP_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0685:'
  write ( *, '(a)' ) '  EXP_VALUES returns values of the exponential function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           EXP(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exp_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test069 ( )

!*****************************************************************************80
!
!! TEST069 tests EXP3_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST069:'
  write ( *, '(a)' ) '  EXP3_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the EXP3 Integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          EXP3_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exp3_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test070 ( )

!*****************************************************************************80
!
!! TEST070 tests EXPONENTIAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST070:'
  write ( *, '(a)' ) '  EXPONENTIAL_CDF_VALUES returns values of'
  write ( *, '(a)' ) '  the Exponential CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          LAMBDA          X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exponential_cdf_values ( n_data, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) lambda, x, fx

  end do

  return
end
subroutine test071 ( )

!*****************************************************************************80
!
!! TEST071 tests EXTREME_VALUES_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST071:'
  write ( *, '(a)' ) '  EXTREME_VALUES_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Extreme Values Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       ALPHA          BETA             X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call extreme_values_cdf_values ( n_data, alpha, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) alpha, beta, x, fx

  end do

  return
end
subroutine test072 ( )

!*****************************************************************************80
!
!! TEST072 tests F_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072:'
  write ( *, '(a)' ) '  F_CDF_VALUES returns values of'
  write ( *, '(a)' ) '  the F cumulative density function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A         B          X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine test073 ( )

!*****************************************************************************80
!
!! TEST073 tests F_NONCENTRAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) fx
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST073:'
  write ( *, '(a)' ) '  F_NONCENTRAL_CDF_VALUES returns values of'
  write ( *, '(a)' ) '  the F cumulative density function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A         B      LAMBDA       X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f10.6,2x,g14.6,2x,g14.6)' ) a, b, lambda, x, fx

  end do

  return
end
subroutine test074 ( )

!*****************************************************************************80
!
!! TEST074 tests FRESNEL_COS_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST074:'
  write ( *, '(a)' ) '  FRESNEL_COS_VALUES returns values of'
  write ( *, '(a)' ) '  the Fresnel cosine integral C(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           C(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call fresnel_cos_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test075 ( )

!*****************************************************************************80
!
!! TEST075 tests FRESNEL_SIN_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075:'
  write ( *, '(a)' ) '  FRESNEL_SIN_VALUES returns values of'
  write ( *, '(a)' ) '  the Fresnel sine integral S(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           S(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call fresnel_sin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0755 ( )

!*****************************************************************************80
!
!! TEST0755 tests FROBENIUS_NUMBER_ORDER2_VALUES.
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
  integer ( kind = 4 ) f
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0755:'
  write ( *, '(a)' ) '  FROBENIUS_NUMBER_ORDER2_VALUES returns values of '
  write ( *, '(a)' ) '  the Frobenius number of order 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        C1        C2     F(C1,C2)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call frobenius_number_order2_values ( n_data, c1, c2, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) c1, c2, f

  end do

  return
end
subroutine test0756 ( )

!*****************************************************************************80
!
!! TEST0756 tests FROBENIUS_NUMBER_ORDER_VALUES, FROBENIUS_NUMBER_DATA_VALUES.
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

  integer ( kind = 4 ), allocatable, dimension ( : ) :: c
  integer ( kind = 4 ) f
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) order

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0756:'
  write ( *, '(a)' ) '  FROBENIUS_NUMBER_ORDER_VALUES returns the order for'
  write ( *, '(a)' ) '  a Frobenius problem;'
  write ( *, '(a)' ) '  FROBENIUS_NUMBER_DATA_VALUES returns the corresponding'
  write ( *, '(a)' ) '  coin denominations.'

  n_data = 0

  do

    call frobenius_number_order_values ( n_data, order )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( c(1:order) )

    call frobenius_number_data_values ( n_data, order, c, f )

    write ( *, '(a)'        ) ' '
    write ( *, '(a,i8)'    ) '  Order = ', order
    write ( *, '(10(2x,i6))' ) c(1:order)
    write ( *, '(a,i12)'   ) '  Frobenius number = ', f

    deallocate ( c )

  end do

  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 tests GAMMA_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076:'
  write ( *, '(a)' ) '  GAMMA_VALUES returns values of the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            GAMMA(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test077 ( )

!*****************************************************************************80
!
!! TEST077 tests GAMMA_CDF_VALUES.
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

  real ( kind = 8 ) fx
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST077:'
  write ( *, '(a)' ) '  GAMMA_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the GAMMA Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       MU             SIGMA            X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 tests GAMMA_INC_VALUES.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078:'
  write ( *, '(a)' ) '  GAMMA_INC_VALUES returns values of '
  write ( *, '(a)' ) '  the incomplete Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               X            GAMMA_INC(A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine test079 ( )

!*****************************************************************************80
!
!! TEST079 tests GAMMA_LOG_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079:'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns values of '
  write ( *, '(a)' ) '  the logarithm of the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            GAMMA_LOG(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test080 ( )

!*****************************************************************************80
!
!! TEST080 tests GEGENBAUER_POLY_VALUES.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST080:'
  write ( *, '(a)' ) '  GEGENBAUER_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Gegenbauer polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       A            X     G(N,A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gegenbauer_poly_values ( n_data, n, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2x,g14.6)' ) n, a, x, fx

  end do

  return
end
subroutine test081 ( )

!*****************************************************************************80
!
!! TEST081 tests GEOMETRIC_CDF_VALUES.
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

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST081:'
  write ( *, '(a)' ) '  GEOMETRIC_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Geometric Probability Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X        P         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call geometric_cdf_values ( n_data, x, p, cdf )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) x, p, cdf

  end do

  return
end
subroutine test082 ( )

!*****************************************************************************80
!
!! TEST082 tests GOODWIN_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST082:'
  write ( *, '(a)' ) '  GOODWIN_VALUES returns values of'
  write ( *, '(a)' ) '  the Goodwin and Staton function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call goodwin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test083 ( )

!*****************************************************************************80
!
!! TEST083 tests GUD_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083:'
  write ( *, '(a)' ) '  GUD_VALUES returns values of '
  write ( *, '(a)' ) '  the Gudermannian function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            GUD(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gud_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0835 ( )

!*****************************************************************************80
!
!! TEST0835 tests HERMITE_FUNCTION_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0835:'
  write ( *, '(a)' ) '  HERMITE_FUNCTION_VALUES returns values of '
  write ( *, '(a)' ) '  the Hermite functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X               Hf(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_function_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test084 ( )

!*****************************************************************************80
!
!! TEST084 tests HERMITE_POLY_PHYS_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST084:'
  write ( *, '(a)' ) '  HERMITE_POLY_PHYS_VALUES returns values of '
  write ( *, '(a)' ) '  the physicist''s Hermite polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X               H(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_poly_phys_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0843 ( )

!*****************************************************************************80
!
!! TEST0843 tests HERMITE_POLY_PROB_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0843:'
  write ( *, '(a)' ) '  HERMITE_POLY_PROB_VALUES returns values of '
  write ( *, '(a)' ) '  the probabilist''s Hermite polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X               He(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_poly_prob_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0845 ( )

!*****************************************************************************80
!
!! TEST0845 tests HYPER_2F1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0845:'
  write ( *, '(a)' ) '  HYPER_2F1_VALUES returns values of '
  write ( *, '(a)' ) '  the hypergeometric 2F1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A           B           C            X       Hyper_2F1(A,B,C,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hyper_2f1_values ( n_data, a, b, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6,2x,g24.16)' ) a, b, c, x, fx

  end do

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests HYPERGEOMETRIC_CDF_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) pop
  integer ( kind = 4 ) sam
  integer ( kind = 4 ) suc
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085:'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Hypergeometric Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       SAM       SUC       POP         X      HyperCDF(S,S,P)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hypergeometric_cdf_values ( n_data, sam, suc, pop, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8,2x,g24.16)' ) sam, suc, pop, x, fx

  end do

  return
end
subroutine test0855 ( )

!*****************************************************************************80
!
!! TEST0855 tests HYPERGEOMETRIC_PDF_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) pop
  integer ( kind = 4 ) sam
  integer ( kind = 4 ) suc
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0855:'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Hypergeometric Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       SAM       SUC       POP         X      HyperPDF(S,S,P)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hypergeometric_pdf_values ( n_data, sam, suc, pop, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8,2x,g24.16)' ) sam, suc, pop, x, fx

  end do

  return
end
subroutine test086 ( )

!*****************************************************************************80
!
!! TEST086 tests FACTORIAL_VALUES.
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

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST086:'
  write ( *, '(a)' ) '  FACTORIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N            N!'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test087 ( )

!*****************************************************************************80
!
!! TEST087 tests FACTORIAL2_VALUES.
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

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST087:'
  write ( *, '(a)' ) '  FACTORIAL2_VALUES returns values of '
  write ( *, '(a)' ) '  the double factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N           N!!'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call factorial2_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test088 ( )

!*****************************************************************************80
!
!! TEST088 tests FACTORIAL_RISING_VALUES.
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

  integer ( kind = 4 ) fmn
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST088:'
  write ( *, '(a)' ) '  FACTORIAL_RISING_VALUES returns some exact values'
  write ( *, '(a)' ) '  of the rising factorial function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M         N      Factorial_rising(M,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call factorial_rising_values ( n_data, m, n, fmn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) m, n, fmn

  end do

  return
end
subroutine test089 ( )

!*****************************************************************************80
!
!! TEST089 tests I0ML0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST089:'
  write ( *, '(a)' ) '  I0ML0_VALUES returns values of'
  write ( *, '(a)' ) '  the I0ML0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i0ml0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test090 ( )

!*****************************************************************************80
!
!! TEST090 tests I1ML1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST090:'
  write ( *, '(a)' ) '  I1ML1_VALUES returns values of'
  write ( *, '(a)' ) '  the I1ML1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i1ml1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test091 ( )

!*****************************************************************************80
!
!! TEST091 tests JACOBI_CN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST091:'
  write ( *, '(a)' ) '  JACOBI_CN_VALUES returns values of '
  write ( *, '(a)' ) '  the Jacobi elliptic CN function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A         X       CN(A,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_cn_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,f10.4,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine test092 ( )

!*****************************************************************************80
!
!! TEST092 tests JACOBI_DN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST092:'
  write ( *, '(a)' ) '  JACOBI_DN_VALUES returns values of '
  write ( *, '(a)' ) '  the Jacobi elliptic DN function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A         X       DN(A,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_dn_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,f10.4,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine test093 ( )

!*****************************************************************************80
!
!! TEST093 tests JACOBI_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST093:'
  write ( *, '(a)' ) '  JACOBI_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Jacobi polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         A         B      X       J(N,A,B)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_poly_values ( n_data, n, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,f24.16,2x,g24.16)' ) n, a, b, x, fx

  end do

  return
end
subroutine test094 ( )

!*****************************************************************************80
!
!! TEST094 tests JACOBI_SN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST094:'
  write ( *, '(a)' ) '  JACOBI_SN_VALUES returns values of '
  write ( *, '(a)' ) '  the Jacobi elliptic SN function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A         X       SN(A,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_sn_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,f10.4,2x,g14.6)' ) a, x, fx

  end do

  return
end
subroutine test0945 ( )

!*****************************************************************************80
!
!! TEST0945 tests JED_CE_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) d
  real ( kind = 8 ) f
  real ( kind = 8 ) jed
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) m
  integer ( kind = 4 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0945:'
  write ( *, '(a)' ) '  JED_CE_VALUES returns:'
  write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
  write ( *, '(a)' ) '  YMDF, the corresponding year, month, day, fraction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED          Y   M   D    F'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_ce_values ( n_data, jed, y, m, d, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,i6,2x,i2,2x,i2,2x,f6.4)' ) jed, y, m, d, f

  end do

  return
end
subroutine test095 ( )

!*****************************************************************************80
!
!! TEST095 tests JED_MJD_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) jed
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) mjd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST095:'
  write ( *, '(a)' ) '  JED_MJD_VALUES returns:'
  write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
  write ( *, '(a)' ) '  MJD, the corresponding Modified Julian Day count.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED           MJD'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_mjd_values ( n_data, jed, mjd )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,f12.2)' ) jed, mjd

  end do

  return
end
subroutine test096 ( )

!*****************************************************************************80
!
!! TEST096 tests JED_RD_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) jed
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) rd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST096:'
  write ( *, '(a)' ) '  JED_RD_VALUES returns:'
  write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
  write ( *, '(a)' ) '  RD, the corresponding Reingold Dershowitz Day count.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED            RD'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_rd_values ( n_data, jed, rd )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,f12.2)' ) jed, rd

  end do

  return
end
subroutine test097 ( )

!*****************************************************************************80
!
!! TEST097 tests JED_WEEKDAY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) jed
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) weekday
  character ( len = 9 ), dimension ( 7 ) :: weekday_name = (/ &
    'Sunday   ', 'Monday   ', 'Tuesday  ', 'Wednesday', 'Thursday ', &
    'Friday   ', 'Saturday ' /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST097:'
  write ( *, '(a)' ) '  JED_WEEKDAY_VALUES returns Julian Ephemeris Dates '
  write ( *, '(a)' ) '  (JED) and the corresponding weekday'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED     #  Weekday'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_weekday_values ( n_data, jed, weekday )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,i1,2x,a9)' ) jed, weekday, weekday_name(weekday)

  end do

  return
end
subroutine test0972 ( )

!*****************************************************************************80
!
!! TEST0972 tests KEI0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0972:'
  write ( *, '(a)' ) '  KEI0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function KEI of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KEI0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call kei0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0973 ( )

!*****************************************************************************80
!
!! TEST0973 tests KEI1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0973:'
  write ( *, '(a)' ) '  KEI1_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function KEI of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KEI1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call kei1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0974 ( )

!*****************************************************************************80
!
!! TEST0974 tests KER0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0974:'
  write ( *, '(a)' ) '  KER0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function KER of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KER0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ker0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test0975 ( )

!*****************************************************************************80
!
!! TEST0975 tests KER1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0975:'
  write ( *, '(a)' ) '  KER1_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function KER of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KER1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ker1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test098 ( )

!*****************************************************************************80
!
!! TEST098 tests LAGUERRE_ASSOCIATED_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST098:'
  write ( *, '(a)' ) '  LAGUERRE_ASSOCIATED_VALUES returns values of '
  write ( *, '(a)' ) '  the associated Laguerre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M      X            L(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laguerre_associated_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine test099 ( )

!*****************************************************************************80
!
!! TEST099 tests LAGUERRE_POLYNOMIAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST099:'
  write ( *, '(a)' ) '  LAGUERRE_POLYNOMIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the Laguerre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            L(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laguerre_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test0995 ( )

!*****************************************************************************80
!
!! TEST0995 tests LAMBERT_W_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0995:'
  write ( *, '(a)' ) '  LAMBERT_W_VALUES returns values of '
  write ( *, '(a)' ) '  the Lambert W function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           W(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lambert_w_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test100 ( )

!*****************************************************************************80
!
!! TEST100 tests LAPLACE_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) beta
  real ( kind = 8 ) fx
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST100:'
  write ( *, '(a)' ) '  LAPLACE_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Laplace CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          MU              BETA            X           CDF(MU,BETA,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laplace_cdf_values ( n_data, mu, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g24.16)' ) mu, beta, x, fx

  end do

  return
end
subroutine test101 ( )

!*****************************************************************************80
!
!! TEST101 tests LEGENDRE_ASSOCIATED_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST101:'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_VALUES returns values of '
  write ( *, '(a)' ) '  the associated Legendre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M          X            P(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine test1015 ( )

!*****************************************************************************80
!
!! TEST1015 tests LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1015:'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES returns values of '
  write ( *, '(a)' ) '  the normalized associated Legendre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M          X            P(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_normalized_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine test1016 ( )

!*****************************************************************************80
!
!! TEST1016 tests LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1016:'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES returns values of '
  write ( *, '(a)' ) '  the associated Legendre polynomials normalized for a sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M          X            P(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_normalized_sphere_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine test102 ( )

!*****************************************************************************80
!
!! TEST102 tests LEGENDRE_POLY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST102:'
  write ( *, '(a)' ) '  LEGENDRE_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Legendre PN polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            P(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine test103 ( )

!*****************************************************************************80
!
!! TEST103 tests LEGENDRE_FUNCTION_Q_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST103:'
  write ( *, '(a)' ) '  LEGENDRE_FUNCTION_Q_VALUES returns values of '
  write ( *, '(a)' ) '  the Legendre QN polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X           Q(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_function_q_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) n, x, fx

  end do

  return
end
subroutine test1035 ( )

!*****************************************************************************80
!
!! TEST1035 tests LERCH_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) s
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1035:'
  write ( *, '(a)' ) '  LERCH_VALUES returns values of '
  write ( *, '(a)' ) '  the Lerch transcendent function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Z            S        A        CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lerch_values ( n_data, z, s, a, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,i8,2x,f10.4,2x,g24.16)' ) z, s, a, fx

  end do

  return
end
subroutine test104 ( )

!*****************************************************************************80
!
!! TEST104 tests LOBACHEVSKY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST104:'
  write ( *, '(a)' ) '  LOBACHEVSKY_VALUES returns values of'
  write ( *, '(a)' ) '  the Lobachevsky function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lobachevsky_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test1037 ( )

!*****************************************************************************80
!
!! TEST1037 tests LOG_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1037:'
  write ( *, '(a)' ) &
    '  LOG_VALUES returns values of the natural logarithm function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            LN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests LOG_NORMAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105:'
  write ( *, '(a)' ) '  LOG_NORMAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Log Normal Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call log_normal_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine test106 ( )

!*****************************************************************************80
!
!! TEST106 tests LOG_SERIES_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST106:'
  write ( *, '(a)' ) '  LOG_SERIES_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Log Series Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     T          N   CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call log_series_cdf_values ( n_data, t, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,i8,2x,g24.16)' ) t, n, fx

  end do

  return
end
subroutine test107 ( )

!*****************************************************************************80
!
!! TEST107 tests LOGARITHMIC_INTEGRAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST107:'
  write ( *, '(a)' ) '  LOGARITHMIC_INTEGAL_VALUES returns values of'
  write ( *, '(a)' ) '  the logarithmic integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            LI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call logarithmic_integral_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test108 ( )

!*****************************************************************************80
!
!! TEST108 tests LOGISTIC_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) beta
  real ( kind = 8 ) fx
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST108:'
  write ( *, '(a)' ) '  LOGISTIC_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Logistic Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      MU     BETA     X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call logistic_cdf_values ( n_data, mu, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, beta, x, fx

  end do

  return
end
subroutine test10825 ( )

!*****************************************************************************80
!
!! TEST10825 tests MATHIEU_EVEN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10825:'
  write ( *, '(a)' ) '  MATHIEU_EVEN_VALUES returns values of '
  write ( *, '(a)' ) '  the eigenvalues of the Mathieu differential'
  write ( *, '(a)' ) '  equation associated with even periodic solutions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R         Q       A(R,Q)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mathieu_even_values ( n_data, r, q, a )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) r, q, a

  end do

  return
end
subroutine test10850 ( )

!*****************************************************************************80
!
!! TEST10850 tests MATHIEU_ODD_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10850:'
  write ( *, '(a)' ) '  MATHIEU_ODD_VALUES returns values of '
  write ( *, '(a)' ) '  the eigenvalues of the Mathieu differential'
  write ( *, '(a)' ) '  equation associated with odd periodic solutions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R         Q       B(R,Q)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mathieu_odd_values ( n_data, r, q, b )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) r, q, b

  end do

  return
end
subroutine test10875 ( )

!*****************************************************************************80
!
!! TEST10875 tests MERTENS_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10875:'
  write ( *, '(a)' ) '  MERTENS_VALUES returns values of '
  write ( *, '(a)' ) '  the Mertens function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N         MERTENS(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mertens_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test109 ( )

!*****************************************************************************80
!
!! TEST109 tests MOEBIUS_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST109:'
  write ( *, '(a)' ) '  MOEBIUS_VALUES returns values of '
  write ( *, '(a)' ) '  the Moebius function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N         MU(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call moebius_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test110 ( )

!*****************************************************************************80
!
!! TEST110 tests NEGATIVE_BINOMIAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) f
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  integer ( kind = 4 ) s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST110:'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Negative Binomial Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       F       S         P         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call negative_binomial_cdf_values ( n_data, f, s, p, cdf )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) f, s, p, cdf

  end do

  return
end
subroutine test1105 ( )

!*****************************************************************************80
!
!! TEST1105 tests NINE_J_VALUES.
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

  real ( kind = 8 ) fx
  real ( kind = 8 ) j1
  real ( kind = 8 ) j2
  real ( kind = 8 ) j3
  real ( kind = 8 ) j4
  real ( kind = 8 ) j5
  real ( kind = 8 ) j6
  real ( kind = 8 ) j7
  real ( kind = 8 ) j8
  real ( kind = 8 ) j9
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1105:'
  write ( *, '(a)' ) '  NINE_J_VALUES returns values of '
  write ( *, '(a)' ) '  the Wigner 9J coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      J4      J5      J6' // &
    '      J7      J8      J9        NINE_J'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call nine_j_values ( n_data, j1, j2, j3, j4, j5, j6, j7, j8, j9, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,' // &
      'f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, j4, j5, j6, j7, j8, j9, fx

  end do

  return
end
subroutine test111 ( )

!*****************************************************************************80
!
!! TEST111 tests NORMAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST111:'
  write ( *, '(a)' ) '  NORMAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Normal Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine test112 ( )

!*****************************************************************************80
!
!! TEST112 tests NORMAL_01_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST112:'
  write ( *, '(a)' ) '  NORMAL_01_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Normal Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test113 ( )

!*****************************************************************************80
!
!! TEST113 tests OMEGA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST113:'
  write ( *, '(a)' ) '  OMEGA_VALUES returns values of '
  write ( *, '(a)' ) '  the Omega function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N           OMEGA(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call omega_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i12,2x,i12)' ) n, fn

  end do

  return
end
subroutine test1135 ( )

!*****************************************************************************80
!
!! TEST1135 tests OWEN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) h
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1135:'
  write ( *, '(a)' ) '  OWEN_VALUES returns values of '
  write ( *, '(a)' ) '  the Owen T function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      H       A       T'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g24.16)' ) h, a, t

  end do

  return
end
subroutine test114 ( )

!*****************************************************************************80
!
!! TEST114 tests PARTITION_COUNT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST114:'
  write ( *, '(a)' ) '  PARTITION_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the integer partition count function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         P(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call partition_count_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test115 ( )

!*****************************************************************************80
!
!! TEST115 tests PARTITION_DISTINCT_COUNT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115:'
  write ( *, '(a)' ) '  PARTITION_DISTINCT_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the integer distinct partition count function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         Q(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call partition_distinct_count_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test116 ( )

!*****************************************************************************80
!
!! TEST116 tests PHI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST116:'
  write ( *, '(a)' ) '  PHI_VALUES returns values of '
  write ( *, '(a)' ) '  the PHI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         PHI(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call phi_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test117 ( )

!*****************************************************************************80
!
!! TEST117 tests PI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST117:'
  write ( *, '(a)' ) '  PI_VALUES returns values of '
  write ( *, '(a)' ) '  the PI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             N         PI(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call pi_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i12,2x,i12)' ) n, fn

  end do

  return
end
subroutine test118 ( )

!*****************************************************************************80
!
!! TEST118 tests POISSON_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST118:'
  write ( *, '(a)' ) '  POISSON_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Poisson Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A          X    CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call poisson_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,i8,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine test1185 ( )

!*****************************************************************************80
!
!! TEST1185 tests POLYLOGARITHM_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1185:'
  write ( *, '(a)' ) '  POLYLOGARITHM_VALUES returns values of '
  write ( *, '(a)' ) '  the polylogarithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      Z        FX'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call polylogarithm_values ( n_data, n, z, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f10.4,2x,g24.16)' ) n, z, fx

  end do

  return
end
subroutine test119 ( )

!*****************************************************************************80
!
!! TEST119 tests PRANDTL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) pr
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST119:'
  write ( *, '(a)' ) '  PRANDTL_VALUES returns values of '
  write ( *, '(a)' ) '  the Prandtl number of water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            P            Pr(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call prandtl_values ( n_data, tc, p, pr )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, pr

  end do

  return
end
subroutine test120 ( )

!*****************************************************************************80
!
!! TEST120 tests PRIME_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) p

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST120:'
  write ( *, '(a)' ) '  PRIME_VALUES returns values of '
  write ( *, '(a)' ) '  the prime function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           N          P[N]'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call prime_values ( n_data, n, p )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i12,2x,i12)' ) n, p

  end do

  return
end
subroutine test121 ( )

!*****************************************************************************80
!
!! TEST121 tests PSAT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) psat
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST121:'
  write ( *, '(a)' ) '  PSAT_VALUES returns values of '
  write ( *, '(a)' ) '  the saturation pressure of water '
  write ( *, '(a)' ) '  as a function of temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            PSAT(T)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psat_values ( n_data, tc, psat )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) tc, psat

  end do

  return
end
subroutine test122 ( )

!*****************************************************************************80
!
!! TEST122 tests PSI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST122:'
  write ( *, '(a)' ) '  PSI_VALUES returns values of '
  write ( *, '(a)' ) '  the PSI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            PSI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test123 ( )

!*****************************************************************************80
!
!! TEST123 tests R8_FACTORIAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST123:'
  write ( *, '(a)' ) '  R8_FACTORIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the factorial function (using real arithmetic).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        N       Factorial(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, fn

  end do

  return
end
subroutine test124 ( )

!*****************************************************************************80
!
!! TEST124 tests R8_FACTORIAL_LOG_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST124:'
  write ( *, '(a)' ) '  R8_FACTORIAL_LOG_VALUES returns values of '
  write ( *, '(a)' ) '  the logarithm of the factorial function '
  write ( *, '(a)' ) '  (using real arithmetic).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        N       Log(Factorial(N))'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_log_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, fn

  end do

  return
end
subroutine test125 ( )

!*****************************************************************************80
!
!! TEST125 tests SECVIR_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) tc
  real ( kind = 8 ) vir

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125:'
  write ( *, '(a)' ) '  SECVIR_VALUES returns values of '
  write ( *, '(a)' ) '  the second virial coefficient of water '
  write ( *, '(a)' ) '  as a function of temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            VIR(T)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call secvir_values ( n_data, tc, vir )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) tc, vir

  end do

  return
end
subroutine test1255 ( )

!*****************************************************************************80
!
!! TEST1255 tests SHI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1255:'
  write ( *, '(a)' ) '  SHI_VALUES returns values of '
  write ( *, '(a)' ) '  the Hyperbolic Sine Integral function SHI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           SHI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call shi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test126 ( )

!*****************************************************************************80
!
!! TEST126 tests SI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST126:'
  write ( *, '(a)' ) '  SI_VALUES returns values of '
  write ( *, '(a)' ) '  the sine integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            SI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call si_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test127 ( )

!*****************************************************************************80
!
!! TEST127 tests SIGMA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST127:'
  write ( *, '(a)' ) '  SIGMA_VALUES returns values of '
  write ( *, '(a)' ) '  the SIGMA function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         SIGMA(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sigma_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test1275 ( )

!*****************************************************************************80
!
!! TEST1275 tests SIN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1275:'
  write ( *, '(a)' ) '  SIN_VALUES returns values of the sine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           SIN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test128 ( )

!*****************************************************************************80
!
!! TEST128 tests SIN_POWER_INT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST128:'
  write ( *, '(a)' ) '  SIN_POWER_INT returns values of '
  write ( *, '(a)' ) '  the integral of SIN(X)^N from A to B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B              N    SIN_POWER_INT(A,B,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sin_power_int_values ( n_data, a, b, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,i8,2x,g14.6)' ) a, b, n, fx

  end do

  return
end
subroutine test1283 ( )

!*****************************************************************************80
!
!! TEST1283 tests SINH_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1283:'
  write ( *, '(a)' ) &
    '  SINH_VALUES returns values of the hyperbolic sine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           SINH(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sinh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test1285 ( )

!*****************************************************************************80
!
!! TEST1285 tests SIX_J_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) j1
  real ( kind = 8 ) j2
  real ( kind = 8 ) j3
  real ( kind = 8 ) j4
  real ( kind = 8 ) j5
  real ( kind = 8 ) j6
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1285:'
  write ( *, '(a)' ) '  SIX_J_VALUES returns values of '
  write ( *, '(a)' ) '  the Wigner 6J coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      J4      J5      J6        SIX_J'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call six_j_values ( n_data, j1, j2, j3, j4, j5, j6, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, j4, j5, j6, fx

  end do

  return
end
subroutine test129 ( )

!*****************************************************************************80
!
!! TEST129 tests SOUND_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST129:'
  write ( *, '(a)' ) '  SOUND_VALUES returns values of '
  write ( *, '(a)' ) '  the spead of sound in water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            P            C(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sound_values ( n_data, tc, p, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, c

  end do

  return
end
subroutine test131 ( )

!*****************************************************************************80
!
!! TEST131 tests SPHERE_UNIT_AREA_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST131:'
  write ( *, '(a)' ) '  SPHERE_UNIT_AREA_VALUES returns values of '
  write ( *, '(a)' ) '  the area of the unit sphere in various dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      N            Area'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sphere_unit_area_values ( n_data, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,g24.16)' ) n, fx

  end do

  return
end
subroutine test132 ( )

!*****************************************************************************80
!
!! TEST132 tests SPHERE_UNIT_VOLUME_VALUES.
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

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST132:'
  write ( *, '(a)' ) '  SPHERE_UNIT_VOLUME_VALUES returns values of '
  write ( *, '(a)' ) '  the volume of the unit sphere in various dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      N            Volume'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sphere_unit_volume_values ( n_data, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,g24.16)' ) n, fx

  end do

  return
end
subroutine test1325 ( )

!*****************************************************************************80
!
!! TEST1325 tests SPHERICAL_HARMONIC_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) phi
  real ( kind = 8 ) theta
  real ( kind = 8 ) yi
  real ( kind = 8 ) yr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1325:'
  write ( *, '(a)' ) '  SPHERICAL_HARMONIC_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical harmonic functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   L   M    THETA       PHI       Yr                         Yi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call spherical_harmonic_values ( n_data, l, m, theta, phi, yr, yi )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i2,2x,i2,2x,f8.4,2x,f8.4,2x,g24.16,2x,g24.16)' ) &
      l, m, theta, phi, yr, yi

  end do

  return
end
subroutine test130 ( )

!*****************************************************************************80
!
!! TEST130 tests SQRT_VALUES.
!
!  Discussion:
!
!    In this example, we suggest how the tabulated values can be
!    used to look for large discrepancies.  (In fact, in this case,
!    we suppose there will be none!).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST130:'
  write ( *, '(a)' ) '  SQRT evaluates the square root function.'
  write ( *, '(a)' ) '  SQRT_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X      Exact F       SQRT(X)        Diff'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sqrt_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = sqrt ( x )

    diff = abs ( fx - fx2 )

    write ( *, '(2x,f14.4,2x,g14.6,2x,g14.6,2x,g14.6)' ) x, fx, fx2, diff

  end do

  return
end
subroutine test133 ( )

!*****************************************************************************80
!
!! TEST133 tests STIRLING1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) m
  integer ( kind = 4 ) s1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST133:'
  write ( *, '(a)' ) '  STIRLING1_VALUES returns values of '
  write ( *, '(a)' ) '  the Stirling numbers of the first kind.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M        S1(N,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call stirling1_values ( n_data, n, m, s1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) n, m, s1

  end do

  return
end
subroutine test134 ( )

!*****************************************************************************80
!
!! TEST134 tests STIRLING2_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) m
  integer ( kind = 4 ) s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST134:'
  write ( *, '(a)' ) '  STIRLING2_VALUES returns values of '
  write ( *, '(a)' ) '  the Stirling numbers of the first kind.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M        S2(N,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call stirling2_values ( n_data, n, m, s2 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) n, m, s2

  end do

  return
end
subroutine test135 ( )

!*****************************************************************************80
!
!! TEST135 tests STROMGEN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST135:'
  write ( *, '(a)' ) '  STROMGEN_VALUES returns values of'
  write ( *, '(a)' ) '  the Stromgen function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call stromgen_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test136 ( )

!*****************************************************************************80
!
!! TEST136 tests STRUVE_H0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST136:'
  write ( *, '(a)' ) '  STRUVE_H0_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve H0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            H0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_h0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test137 ( )

!*****************************************************************************80
!
!! TEST137 tests STRUVE_H1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST137:'
  write ( *, '(a)' ) '  STRUVE_H1_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve H1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            H1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_h1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test138 ( )

!*****************************************************************************80
!
!! TEST138 tests STRUVE_L0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST138:'
  write ( *, '(a)' ) '  STRUVE_L0_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve L0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            L0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_l0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test139 ( )

!*****************************************************************************80
!
!! TEST139 tests STRUVE_L1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST139:'
  write ( *, '(a)' ) '  STRUVE_L1_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve L1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            L1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_l1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test140 ( )

!*****************************************************************************80
!
!! TEST140 tests STUDENT_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 November 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) c
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST140:'
  write ( *, '(a)' ) '  STUDENT_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Student T Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      C     X       CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_cdf_values ( n_data, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,f10.4,2x,g24.16)' ) c, x, fx

  end do

  return
end
subroutine test141 ( )

!*****************************************************************************80
!
!! TEST141 tests STUDENT_NONCENTRAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) df
  real ( kind = 8 ) fx
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST141:'
  write ( *, '(a)' ) '  STUDENT_NONCENTRAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the noncentral Student T Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      LAMBDA       DF     CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.4,2x,f10.4,2x,i8,2x,g14.6)' ) x, lambda, df, fx

  end do

  return
end
subroutine test1415 ( )

!*****************************************************************************80
!
!! TEST1415 tests SUBFACTORIAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1415:'
  write ( *, '(a)' ) '  SUBFACTORIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the subfactorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N     Subfactorial[N]'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call subfactorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test142 ( )

!*****************************************************************************80
!
!! TEST142 tests SURTEN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST142:'
  write ( *, '(a)' ) '  SURTEN_VALUES returns values of '
  write ( *, '(a)' ) '  the surface tension of water '
  write ( *, '(a)' ) '  as a function of temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            SIGMA(T)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call surten_values ( n_data, tc, sigma )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) tc, sigma

  end do

  return
end
subroutine test143 ( )

!*****************************************************************************80
!
!! TEST143 tests SYNCH1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST143:'
  write ( *, '(a)' ) '  SYNCH1_VALUES returns values of '
  write ( *, '(a)' ) '  the synchrotron radiation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            S1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call synch1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test144 ( )

!*****************************************************************************80
!
!! TEST144 tests SYNCH2_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST144:'
  write ( *, '(a)' ) '  SYNCH2_VALUES returns values of '
  write ( *, '(a)' ) '  the synchrotron radiation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            S2(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call synch2_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test1445 ( )

!*****************************************************************************80
!
!! TEST1445 tests TAN_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1445:'
  write ( *, '(a)' ) '  TAN_VALUES returns values of the tangent function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           TAN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tan_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test1447 ( )

!*****************************************************************************80
!
!! TEST1447 tests TANH_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1447:'
  write ( *, '(a)' ) &
    '  TANH_VALUES returns values of the hyperbolic tangent function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           TANH(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tanh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test145 ( )

!*****************************************************************************80
!
!! TEST145 tests TAU_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) fn
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST145:'
  write ( *, '(a)' ) '  TAU_VALUES returns values of '
  write ( *, '(a)' ) '  the TAU function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         TAU(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tau_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine test146 ( )

!*****************************************************************************80
!
!! TEST146 tests THERCON_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST146:'
  write ( *, '(a)' ) '  THERCON_VALUES returns values of '
  write ( *, '(a)' ) '  the thermal conductivity of water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            P            LAMBDA(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call thercon_values ( n_data, tc, p, lambda )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, lambda

  end do

  return
end
subroutine test1465 ( )

!*****************************************************************************80
!
!! TEST1465 tests THREE_J_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) j1
  real ( kind = 8 ) j2
  real ( kind = 8 ) j3
  real ( kind = 8 ) m1
  real ( kind = 8 ) m2
  real ( kind = 8 ) m3
  integer ( kind = 4 ) n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1465:'
  write ( *, '(a)' ) '  THREE_J_VALUES returns values of '
  write ( *, '(a)' ) '  the Wigner 3J coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      M1      M2      M3        THREE_J'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call three_j_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, m1, m2, m3, fx

  end do

  return
end
subroutine test147 ( )

!*****************************************************************************80
!
!! TEST147 tests TRAN02_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST147:'
  write ( *, '(a)' ) '  TRAN02_VALUES returns values of '
  write ( *, '(a)' ) '  the order 2 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T2(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran02_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test148 ( )

!*****************************************************************************80
!
!! TEST148 tests TRAN03_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST148:'
  write ( *, '(a)' ) '  TRAN03_VALUES returns values of '
  write ( *, '(a)' ) '  the order 3 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T3(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran03_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test149 ( )

!*****************************************************************************80
!
!! TEST149 tests TRAN04_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST149:'
  write ( *, '(a)' ) '  TRAN04_VALUES returns values of '
  write ( *, '(a)' ) '  the order 4 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T4(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran04_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test150 ( )

!*****************************************************************************80
!
!! TEST150 tests TRAN05_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST150:'
  write ( *, '(a)' ) '  TRAN05_VALUES returns values of '
  write ( *, '(a)' ) '  the order 5 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T5(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran05_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test151 ( )

!*****************************************************************************80
!
!! TEST151 tests TRAN06_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST151:'
  write ( *, '(a)' ) '  TRAN06_VALUES returns values of '
  write ( *, '(a)' ) '  the order 6 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T6(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran06_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test152 ( )

!*****************************************************************************80
!
!! TEST152 tests TRAN07_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST152:'
  write ( *, '(a)' ) '  TRAN07_VALUES returns values of '
  write ( *, '(a)' ) '  the order 7 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T7(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran07_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test153 ( )

!*****************************************************************************80
!
!! TEST153 tests TRAN08_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST153:'
  write ( *, '(a)' ) '  TRAN08_VALUES returns values of '
  write ( *, '(a)' ) '  the order 8 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T8(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran08_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test154 ( )

!*****************************************************************************80
!
!! TEST154 tests TRAN09_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST154:'
  write ( *, '(a)' ) '  TRAN09_VALUES returns values of '
  write ( *, '(a)' ) '  the order 9 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T9(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran09_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine test1545 ( )

!*****************************************************************************80
!
!! TEST1545 tests TRIGAMMA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1545:'
  write ( *, '(a)' ) '  TRIGAMMA_VALUES returns values of '
  write ( *, '(a)' ) '  the TriGamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call trigamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine test155 ( )

!*****************************************************************************80
!
!! TEST155 tests TSAT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST155:'
  write ( *, '(a)' ) '  TSAT_VALUES returns values of '
  write ( *, '(a)' ) '  the saturation temperature '
  write ( *, '(a)' ) '  as a function of pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      P           Tsat(P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tsat_values ( n_data, p, tc )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,g14.6)' ) p, tc

  end do

  return
end
subroutine test156 ( )

!*****************************************************************************80
!
!! TEST156 tests VAN_DER_CORPUT_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) base
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) seed
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST156:'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_VALUES returns values of '
  write ( *, '(a)' ) '  the van der Corput sequence in a given base.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      BASE      SEED    VDC(BASE,SEED)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call van_der_corput_values ( n_data, base, seed, value )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g16.8)' ) base, seed, value

  end do

  return
end
subroutine test157 ( )

!*****************************************************************************80
!
!! TEST157 tests VISCOSITY_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) eta
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) p
  real ( kind = 8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST157:'
  write ( *, '(a)' ) '  VISCOSITY_VALUES returns values of '
  write ( *, '(a)' ) '  the viscosity of water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          T               P          ETA(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call viscosity_values ( n_data, tc, p, eta )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,g14.6)' ) tc, p, eta

  end do

  return
end
subroutine test1575 ( )

!*****************************************************************************80
!
!! TEST1575 tests VON_MISES_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1575:'
  write ( *, '(a)' ) '  VON_MISES_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the von Mises CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B               X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call von_mises_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.6,2x,f14.6,2x,f14.6,2x,g14.6)' ) a, b, x, fx

  end do

  return
end
subroutine test158 ( )

!*****************************************************************************80
!
!! TEST158 tests WEIBULL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) fx
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST158:'
  write ( *, '(a)' ) '  WEIBULL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Weibull Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       ALPHA          BETA             X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call weibull_cdf_values ( n_data, alpha, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) alpha, beta, x, fx

  end do

  return
end
subroutine test159 ( )

!*****************************************************************************80
!
!! TEST159 tests ZETA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) zeta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST159:'
  write ( *, '(a)' ) '  ZETA_VALUES returns values of '
  write ( *, '(a)' ) '  the Riemann Zeta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       ZETA(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call zeta_values ( n_data, n, zeta )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, zeta

  end do

  return
end
