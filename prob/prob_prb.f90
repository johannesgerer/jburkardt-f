program main

!*****************************************************************************80
!
!! MAIN is the main program for PROB_PRB.
!
!  Discussion:
!
!    PROB_PRB calls sample problems for the PROB routines.
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
  write ( *, '(a)' ) 'PROB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the PROB library.'

  call test001 ( )
  call test002 ( )
  call test003 ( )
  call test004 ( )
  call test005 ( )
  call test006 ( )
  call test007 ( )
  call test008 ( )
  call test009 ( )

  call test010 ( )
  call test0105 ( )
  call test0106 ( )
  call test011 ( )
  call test012 ( )
  call test013 ( )
  call test014 ( )
  call test015 ( )
  call test016 ( )

  call test020 ( )
  call test021 ( )
  call test022 ( )
  call test023 ( )
  call test0235 ( )
  call test024 ( )
  call test025 ( )
  call test0251 ( )
  call test0252 ( )
  call test0253 ( )
  call test0254 ( )
  call test026 ( )
  call test027 ( )
  call test0275 ( )
  call test0276 ( )
  call test028 ( )
  call test029 ( )

  call test030 ( )
  call test031 ( )
  call test032 ( )
  call test033 ( )
  call test034 ( )
  call test035 ( )
  call test036 ( )
  call test037 ( )
  call test0375 ( )
  call test038 ( )
  call test039 ( )
  call test0395 ( )

  call test040 ( )
  call test041 ( )
  call test042 ( )
  call test043 ( )
  call test044 ( )
  call test045 ( )
  call test046 ( )
  call test047 ( )
  call test048 ( )
  call test049 ( )

  call test050 ( )
  call test051 ( )
  call test052 ( )
  call test053 ( )
  call test054 ( )
  call test055 ( )
  call test056 ( )
  call test0563 ( )
  call test0564 ( )
  call test0565 ( )
  call test0566 ( )
  call test057 ( )
  call test058 ( )
  call test059 ( )

  call test060 ( )
  call test061 ( )
  call test062 ( )
  call test063 ( )
  call test064 ( )
  call test065 ( )
  call test066 ( )
  call test067 ( )
  call test068 ( )
  call test069 ( )

  call test070 ( )
  call test07025 ( )
  call test0705 ( )
  call test071 ( )
  call test072 ( )
  call test073 ( )
  call test074 ( )
  call test0744 ( )
  call test0745 ( )
  call test075 ( )
  call test076 ( )
  call test077 ( )
  call test078 ( )
  call test079 ( )

  call test080 ( )
  call test081 ( )
  call test082 ( )
  call test083 ( )
  call test084 ( )
  call test085 ( )
  call test086 ( )
  call test087 ( )
  call test088 ( )
  call test089 ( )

  call test090 ( )
  call test091 ( )
  call test092 ( )
  call test093 ( )
  call test094 ( )
  call test095 ( )
  call test096 ( )
  call test0965 ( )
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
  call test109 ( )

  call test110 ( )
  call test111 ( )
  call test112 ( )
  call test113 ( )
  call test114 ( )
  call test1145 ( )
  call test1146 ( )
  call test115 ( )
  call test116 ( )
  call test117 ( )
  call test118 ( )
  call test119 ( )

  call test120 ( )
  call test123 ( )
  call test124 ( )
  call test125 ( )
  call test126 ( )
  call test127 ( )
  call test128 ( )
  call test129 ( )

  call test130 ( )
  call test1304 ( )
  call test1306 ( )
  call test131 ( )
  call test132 ( )
  call test133 ( )
  call test134 ( )
  call test1341 ( )
  call test1342 ( )
  call test1344 ( )
  call test135 ( )
  call test136 ( )
  call test137 ( )
  call test138 ( )
  call test139 ( )

  call test140 ( )
  call test141 ( )
  call test142 ( )
  call test1425 ( )
  call test143 ( )
  call test144 ( )
  call test145 ( )
  call test146 ( )
  call test147 ( )
  call test148 ( )
  call test1485 ( )
  call test1486 ( )
  call test149 ( )

  call test150 ( )
  call test151 ( )
  call test152 ( )
  call test153 ( )
  call test154 ( )
  call test155 ( )
  call test1555 ( )
  call test156 ( )
  call test157 ( )
  call test158 ( )
  call test159 ( )

  call test160 ( )
  call test161 ( )
  call test162 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PROB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests ANGLE_CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) n
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  For the ANGLE PDF:'
  write ( *, '(a)' ) '  ANGLE_CDF evaluates the CDF;'

  n = 5
  x = 0.50D+00

  call angle_cdf ( x, n, cdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter N = ', n
  write ( *, '(a,g14.6)' ) '  PDF argument X =  ', x
  write ( *, '(a,g14.6)' ) '  CDF value =       ', cdf

  return
end
subroutine test002 ( )

!*****************************************************************************80
!
!! TEST002 tests ANGLE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST002'
  write ( *, '(a)' ) '  For the ANGLE PDF:'
  write ( *, '(a)' ) '  ANGLE_PDF evaluates the PDF;'

  n = 5
  x = 0.50D+00

  call angle_pdf ( x, n, pdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter N = ', n
  write ( *, '(a,g14.6)' ) '  PDF argument X =  ', x
  write ( *, '(a,g14.6)' ) '  PDF value =       ', pdf

  return
end
subroutine test003 ( )

!*****************************************************************************80
!
!! TEST003 tests ANGLE_MEAN;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) mean
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST003'
  write ( *, '(a)' ) '  For the ANGLE PDF:'
  write ( *, '(a)' ) '  ANGLE_MEAN computes the mean;'

  n = 5
  call angle_mean ( n, mean )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter N = ', n
  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean

  return
end
subroutine test004 ( )

!*****************************************************************************80
!
!! TEST004 tests ANGLIT_CDF, ANGLIT_CDF_INV, ANGLIT_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST004'
  write ( *, '(a)' ) '  For the Anglit PDF:'
  write ( *, '(a)' ) '  ANGLIT_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  ANGLIT_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  ANGLIT_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call anglit_sample ( seed, x )

    call anglit_pdf ( x, pdf )

    call anglit_cdf ( x, cdf )

    call anglit_cdf_inv ( cdf, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests ANGLIT_MEAN, ANGLIT_SAMPLE, ANGLIT_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  For the Anglit PDF:'
  write ( *, '(a)' ) '  ANGLIT_MEAN computes the mean;'
  write ( *, '(a)' ) '  ANGLIT_SAMPLE samples;'
  write ( *, '(a)' ) '  ANGLIT_VARIANCE computes the variance.'

  call anglit_mean ( mean )
  call anglit_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call anglit_sample ( seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test006 ( )

!*****************************************************************************80
!
!! TEST006 tests ARCSIN_CDF, ARCSIN_CDF_INV, ARCSIN_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  logical arcsin_check
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST006'
  write ( *, '(a)' ) '  For the Arcsin PDF:'
  write ( *, '(a)' ) '  ARCSIN_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  ARCSIN_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  ARCSIN_PDF evaluates the PDF;'

  a = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a

  if ( .not. arcsin_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call arcsin_sample ( a, seed, x )

    call arcsin_pdf ( x, a, pdf )

    call arcsin_cdf ( x, a, cdf )

    call arcsin_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test007 ( )

!*****************************************************************************80
!
!! TEST007 tests ARCSIN_MEAN, ARCSIN_SAMPLE, ARCSIN_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  logical arcsin_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST007'
  write ( *, '(a)' ) '  For the Arcsin PDF:'
  write ( *, '(a)' ) '  ARCSIN_MEAN computes the mean;'
  write ( *, '(a)' ) '  ARCSIN_SAMPLE samples;'
  write ( *, '(a)' ) '  ARCSIN_VARIANCE computes the variance.'

  do i = 1, 2

    if ( i == 1 ) then
      a = 1.0D+00
    else if ( i == 2 ) then
      a = 16.0D+00
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a

    if ( .not. arcsin_check ( a ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Fatal error!'
      write ( *, '(a)' ) '  The parameters are not legal.'
      return
    end if

    call arcsin_mean ( a, mean )
    call arcsin_variance ( a, variance )

    write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
    write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

    do j = 1, sample_num
      call arcsin_sample ( a, seed, x(j) )
    end do

    call r8vec_mean ( sample_num, x, mean )
    call r8vec_variance ( sample_num, x, variance )
    call r8vec_max ( sample_num, x, imax, xmax )
    call r8vec_min ( sample_num, x, imin, xmin )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
    write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
    write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
    write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
    write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  end do

  return
end
subroutine test008 ( )

!*****************************************************************************80
!
!! TEST008 tests BENFORD_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST008'
  write ( *, '(a)' ) '  For the Benford PDF:'
  write ( *, '(a)' ) '  BENFORD_PDF evaluates the PDF.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N    PDF(N)'
  write ( *, '(a)' ) ' '

  do n = 1, 19

    call benford_pdf ( n, pdf )
    write ( *, '(i8,g14.6)' ) n, pdf

  end do

  return
end
subroutine test009 ( )

!*****************************************************************************80
!
!! TEST009 tests BERNOULLI_CDF, BERNOULLI_CDF_INV, BERNOULLI_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  logical bernoulli_check
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST009'
  write ( *, '(a)' ) '  For the Bernoulli PDF,'
  write ( *, '(a)' ) '  BERNOULLI_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  BERNOULLI_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  BERNOULLI_PDF evaluates the PDF;'

  a = 0.75D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a

  if ( .not. bernoulli_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call bernoulli_sample ( a, seed, x )

    call bernoulli_pdf ( x, a, pdf )

    call bernoulli_cdf ( x, a, cdf )

    call bernoulli_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test010 ( )

!*****************************************************************************80
!
!! TEST010 tests BERNOULLI_MEAN, BERNOULLI_SAMPLE, BERNOULLI_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  logical bernoulli_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST010'
  write ( *, '(a)' ) '  For the Bernoulli PDF:'
  write ( *, '(a)' ) '  BERNOULLI_MEAN computes the mean;'
  write ( *, '(a)' ) '  BERNOULLI_SAMPLE samples;'
  write ( *, '(a)' ) '  BERNOULLI_VARIANCE computes the variance.'

  a = 0.75D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a

  if ( .not. bernoulli_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call bernoulli_mean ( a, mean )
  call bernoulli_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call bernoulli_sample ( a, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test0105 ( )

!*****************************************************************************80
!
!! TEST0105 demonstrates the use of BESSEL_I0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bessel_i0
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0105:'
  write ( *, '(a)' ) '  BESSEL_I0 computes values of '
  write ( *, '(a)' ) '    the Bessel I0 function.'
  write ( *, '(a)' ) '  BESSEL_I0_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X            Exact                  BESSEL_I0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = bessel_i0 ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine test0106 ( )

!*****************************************************************************80
!
!! TEST0106 demonstrates the use of BESSEL_I1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bessel_i1
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0106:'
  write ( *, '(a)' ) '  BESSEL_I1 computes values of '
  write ( *, '(a)' ) '    the Bessel I1 function.'
  write ( *, '(a)' ) '  BESSEL_I1_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X            Exact                  BESSEL_I1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = bessel_i1 ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine test011 ( )

!*****************************************************************************80
!
!! TEST011 tests BETA and GAMMA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta
  real ( kind = 8 ) beta1
  real ( kind = 8 ) beta2
  real ( kind = 8 ) gamma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST011'
  write ( *, '(a)' ) '  BETA evaluates the Beta function;'
  write ( *, '(a)' ) '  GAMMA evaluates the Gamma function.'

  a = 2.2D+00
  b = 3.7D+00

  beta1 = beta ( a, b )
  beta2 = gamma ( a ) * gamma ( b ) / gamma ( a + b )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Argument A =                   ', a
  write ( *, '(a,g14.6)' ) '  Argument B =                   ', b
  write ( *, '(a,g14.6)' ) '  Beta(A,B) =                    ', beta1
  write ( *, '(a)' ) '  (Expected value = 0.0454 )'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Gamma(A)*Gamma(B)/Gamma(A+B) = ', beta2

  return
end
subroutine test012 ( )

!*****************************************************************************80
!
!! TEST012 tests BETA_CDF, BETA_CDF_INV and BETA_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical beta_check
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST012'
  write ( *, '(a)' ) '  For the Beta PDF:'
  write ( *, '(a)' ) '  BETA_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  BETA_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  BETA_PDF evaluates the PDF;'

  a = 12.0D+00
  b = 12.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. beta_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call beta_sample ( a, b, seed, x )

    call beta_pdf ( x, a, b, pdf )

    call beta_cdf ( x, a, b, cdf )

    call beta_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test013 ( )

!*****************************************************************************80
!
!! TEST013 tests BETA_INC and BETA_INC_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) beta_inc
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST013:'
  write ( *, '(a)' ) '  BETA_INC evaluates the normalized incomplete Beta'
  write ( *, '(a)' ) '    function BETA_INC(A,B,X).'
  write ( *, '(a)' ) '  BETA_INC_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A      B       X       Exact F       BETA_INC(A,B,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = beta_inc ( a, b, x )

    write ( *, '(2x,3f8.4,2g14.6)' ) a, b, x, fx, fx2

  end do

  return
end
subroutine test014 ( )

!*****************************************************************************80
!
!! TEST014 tests BETA_MEAN, BETA_SAMPLE and BETA_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical beta_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST014'
  write ( *, '(a)' ) '  For the Beta PDF:'
  write ( *, '(a)' ) '  BETA_MEAN computes the mean;'
  write ( *, '(a)' ) '  BETA_SAMPLE samples;'
  write ( *, '(a)' ) '  BETA_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. beta_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call beta_mean ( a, b, mean )
  call beta_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call beta_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test015 ( )

!*****************************************************************************80
!
!! TEST015 tests BETA_BINOMIAL_CDF, BETA_BINOMIAL_CDF_INV and BETA_BINOMIAL_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical beta_binomial_check
  integer ( kind = 4 ) c
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST015'
  write ( *, '(a)' ) '  For the Beta Binomial PDF,'
  write ( *, '(a)' ) '  BETA_BINOMIAL_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  BETA_BINOMIAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  BETA_BINOMIAL_PDF evaluates the PDF;'

  a = 2.0D+00
  b = 3.0D+00
  c = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,i8)'    ) '  PDF parameter C = ', c

  if ( .not. beta_binomial_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call beta_binomial_sample ( a, b, c, seed, x )

    call beta_binomial_pdf ( x, a, b, c, pdf )

    call beta_binomial_cdf ( x, a, b, c, cdf )

    call beta_binomial_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test016 ( )

!*****************************************************************************80
!
!! TEST016 tests BETA_BINOMIAL_MEAN, BETA_BINOMIAL_SAMPLE, BETA_BINOMIAL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical beta_binomial_check
  integer ( kind = 4 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST016'
  write ( *, '(a)' ) '  For the Beta Binomial PDF:'
  write ( *, '(a)' ) '  BETA_BINOMIAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  BETA_BINOMIAL_SAMPLE samples;'
  write ( *, '(a)' ) '  BETA_BINOMIAL_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00
  c = 4

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,i8)'    ) '  PDF parameter C = ', c

  if ( .not. beta_binomial_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call beta_binomial_mean ( a, b, c, mean )
  call beta_binomial_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call beta_binomial_sample ( a, b, c, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test020 ( )

!*****************************************************************************80
!
!! TEST020 tests BINOMIAL_CDF and BINOMIAL_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST020:'
  write ( *, '(a)' ) '  BINOMIAL_CDF evaluates the cumulative distribution'
  write ( *, '(a)' ) '    function for the discrete binomial probability'
  write ( *, '(a)' ) '    density function.'
  write ( *, '(a)' ) '  BINOMIAL_CDF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A is the number of trials;'
  write ( *, '(a)' ) '  B is the probability of success on one trial;'
  write ( *, '(a)' ) '  X is the number of successes;'
  write ( *, '(a)' ) '  BINOMIAL_CDF is the probability of having up to X'
  write ( *, '(a)' ) '  successes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A     B         X   Exact F     BINOMIAL_CDF(A,B,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call binomial_cdf ( x, a, b, fx2 )

    write ( *, '(2x,i8,2x,f8.4,2x,i8,g14.6,g14.6)' ) a, b, x, fx, fx2

  end do

  return
end
subroutine test021 ( )

!*****************************************************************************80
!
!! TEST021 tests BINOMIAL_CDF, BINOMIAL_CDF_INV, BINOMIAL_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  logical binomial_check
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST021'
  write ( *, '(a)' ) '  For the Binomial PDF:'
  write ( *, '(a)' ) '  BINOMIAL_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  BINOMIAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  BINOMIAL_PDF evaluates the PDF;'

  a = 5
  b = 0.65D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. binomial_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call binomial_sample ( a, b, seed, x )

    call binomial_pdf ( x, a, b, pdf )

    call binomial_cdf ( x, a, b, cdf )

    call binomial_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test022 ( )

!*****************************************************************************80
!
!! TEST022 tests BINOMIAL_COEF, BINOMIAL_COEF_LOG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) cnk1
  real ( kind = 8 ) cnk2_log
  real ( kind = 8 ) cnk2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST022'
  write ( *, '(a)' ) '  BINOMIAL_COEF evaluates binomial coefficients.'
  write ( *, '(a)' ) '  BINOMIAL_COEF_LOG evaluates the logarithm.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N     K       C(N,K)'
  write ( *, '(a)' ) ' '

  do n = 0, 4
    do k = 0, n
      call binomial_coef ( n, k, cnk1 )
      call binomial_coef_log ( n, k, cnk2_log )
      cnk2 = exp ( cnk2_log )
      write ( *, '(3i8,g14.6)' ) n, k, cnk1, cnk2
    end do
  end do

  return
end
subroutine test023 ( )

!*****************************************************************************80
!
!! TEST023 tests BINOMIAL_MEAN, BINOMIAL_SAMPLE, BINOMIAL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  logical binomial_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST023'
  write ( *, '(a)' ) '  For the Binomial PDF:'
  write ( *, '(a)' ) '  BINOMIAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  BINOMIAL_SAMPLE samples;'
  write ( *, '(a)' ) '  BINOMIAL_VARIANCE computes the variance.'

  a = 5
  b = 0.30D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. binomial_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call binomial_mean ( a, b, mean )
  call binomial_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call binomial_sample ( a, b, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test0235 ( )

!*****************************************************************************80
!
!! TEST0235 tests BIRTHDAY_CDF, BIRTHDAY_CDF_INV, BIRTHDAY_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n2
  real ( kind = 8 ) pdf

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0235'
  write ( *, '(a)' ) '  For the Birthday PDF,'
  write ( *, '(a)' ) '  BIRTHDAY_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  BIRTHDAY_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  BIRTHDAY_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do n = 1, 30

    call birthday_pdf ( n, pdf )

    call birthday_cdf ( n, cdf )

    call birthday_cdf_inv ( cdf, n2 )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,i8)' ) n, pdf, cdf, n2

  end do

  return
end
subroutine test024 ( )

!*****************************************************************************80
!
!! TEST024 tests BRADFORD_CDF, BRADFORD_CDF_INV, BRADFORD_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical bradford_check
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST024'
  write ( *, '(a)' ) '  For the Bradford PDF:'
  write ( *, '(a)' ) '  BRADFORD_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  BRADFORD_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  BRADFORD_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C = ', c

  if ( .not. bradford_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call bradford_sample ( a, b, c, seed, x )

    call bradford_pdf ( x, a, b, c, pdf )

    call bradford_cdf ( x, a, b, c, cdf )

    call bradford_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test025 ( )

!*****************************************************************************80
!
!! TEST025 tests BRADFORD_MEAN, BRADFORD_SAMPLE, BRADFORD_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical bradford_check
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST025'
  write ( *, '(a)' ) '  For the Bradford PDF:'
  write ( *, '(a)' ) '  BRADFORD_MEAN computes the mean;'
  write ( *, '(a)' ) '  BRADFORD_SAMPLE samples;'
  write ( *, '(a)' ) '  BRADFORD_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C = ', c

  if ( .not. bradford_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call bradford_mean ( a, b, c, mean )
  call bradford_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call bradford_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test0251 ( )

!*****************************************************************************80
!
!! TEST0251 tests BUFFON_LAPLACE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) l
  real ( kind = 8 ) pdf

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0251'
  write ( *, '(a)' ) '  BUFFON_LAPLACE_PDF evaluates the Buffon-Laplace PDF,'
  write ( *, '(a)' ) '  the probability that, on a grid of cells of width A'
  write ( *, '(a)' ) '  and height B, a needle of length L, dropped at random,'
  write ( *, '(a)' ) '  will cross at least one grid line.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A         B         L        PDF'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    a = real ( i, kind = 8 )
    do j = 1, 5
      b = real ( j, kind = 8 )

      do k = 0, 5
        l = real ( k, kind = 8 ) * min ( a, b ) / 5.0D+00
        call buffon_laplace_pdf ( a, b, l, pdf )
        write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) a, b, l, pdf
      end do

      write ( *, '(a)' ) ' '

    end do

  end do

  return
end
subroutine test0252 ( )

!*****************************************************************************80
!
!! TEST0252 tests BUFFON_LAPLACE_SIMULATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) buffon_laplace_simulate
  real ( kind = 8 ) err
  integer ( kind = 4 ) hits
  real ( kind = 8 ) l
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00
  real ( kind = 8 ) pi_est
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) trial_num
  integer ( kind = 4 ), dimension ( test_num ) :: trial_num_test = (/ &
    10, 100, 10000, 1000000 /)

  a = 1.0D+00
  b = 1.0D+00
  l = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0252'
  write ( *, '(a)' ) '  BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace'
  write ( *, '(a)' ) '  needle dropping experiment.  On a grid of cells of '
  write ( *, '(a)' ) '  width A and height B, a needle of length L is dropped'
  write ( *, '(a)' ) '  at random.  We count the number of times it crosses'
  write ( *, '(a)' ) '  at least one grid line, and use this to estimate '
  write ( *, '(a)' ) '  the value of PI.'

  seed = 123456789

  call random_initialize ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  Cell width A =    ', a
  write ( *, '(a,f14.6)' ) '  Cell height B =   ', b
  write ( *, '(a,f14.6)' ) '  Needle length L = ', l
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Trials      Hits          Est(Pi)     Err'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    trial_num = trial_num_test(test)

    hits = buffon_laplace_simulate ( a, b, l, trial_num )

    if ( 0 < hits ) then
      pi_est = ( 2.0D+00 * l * ( a + b ) - l * l ) &
        * real ( trial_num, kind = 8 ) &
        / ( a * b * real ( hits, kind = 8 ) )
    else
      pi_est = huge ( pi_est )
    end if

    err = abs ( pi_est - pi )

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g14.6)' ) trial_num, hits, pi_est, err

  end do

  return
end
subroutine test0253 ( )

!*****************************************************************************80
!
!! TEST0253 tests BUFFON_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) l
  real ( kind = 8 ) pdf

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0253'
  write ( *, '(a)' ) '  BUFFON_PDF evaluates the Buffon PDF,'
  write ( *, '(a)' ) '  the probability that, on a grid of cells of width A,'
  write ( *, '(a)' ) '  a needle of length L, dropped at random,'
  write ( *, '(a)' ) '  will cross at least one grid line.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A         L        PDF'
  write ( *, '(a)' ) ' '

  do i = 1, 5
    a = real ( i, kind = 8 )

    do k = 0, 5
      l = real ( k, kind = 8 ) * a / 5.0D+00
      call buffon_pdf ( a, l, pdf )
      write ( *, '(2x,f8.4,2x,f8.4,2x,g14.6)' ) a, l, pdf
    end do

    write ( *, '(a)' ) ' '

  end do

  return
end
subroutine test0254 ( )

!*****************************************************************************80
!
!! TEST0254 tests BUFFON_SIMULATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 4

  real ( kind = 8 ) a
  integer ( kind = 4 ) buffon_simulate
  real ( kind = 8 ) err
  integer ( kind = 4 ) hits
  real ( kind = 8 ) l
  real ( kind = 8 ), parameter :: pi = 3.141592653589793238462643D+00
  real ( kind = 8 ) pi_est
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  integer ( kind = 4 ) trial_num
  integer ( kind = 4 ), dimension ( test_num ) :: trial_num_test = (/ &
    10, 100, 10000, 1000000 /)

  a = 1.0D+00
  l = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0254'
  write ( *, '(a)' ) '  BUFFON_SIMULATE simulates a Buffon-Laplace'
  write ( *, '(a)' ) '  needle dropping experiment.  On a grid of cells of '
  write ( *, '(a)' ) '  width A, a needle of length L is dropped'
  write ( *, '(a)' ) '  at random.  We count the number of times it crosses'
  write ( *, '(a)' ) '  at least one grid line, and use this to estimate '
  write ( *, '(a)' ) '  the value of PI.'

  seed = 123456789

  call random_initialize ( seed )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.6)' ) '  Cell width A =    ', a
  write ( *, '(a,f14.6)' ) '  Needle length L = ', l
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Trials      Hits          Est(Pi)     Err'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    trial_num = trial_num_test(test)

    hits = buffon_simulate ( a, l, trial_num )

    if ( 0 < hits ) then
      pi_est = ( 2.0D+00 * l ) * real ( trial_num, kind = 8 ) &
        / ( a * real ( hits, kind = 8 ) )
    else
      pi_est = huge ( pi_est )
    end if

    err = abs ( pi_est - pi )

    write ( *, '(2x,i8,2x,i8,2x,f14.6,2x,g14.6)' ) trial_num, hits, pi_est, err

  end do

  return
end
subroutine test026 ( )

!*****************************************************************************80
!
!! TEST026 tests BURR_CDF, BURR_CDF_INV, BURR_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical burr_check
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST026'
  write ( *, '(a)' ) '  For the Burr PDF:'
  write ( *, '(a)' ) '  BURR_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  BURR_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  BURR_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00
  d = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C = ', c
  write ( *, '(a,g14.6)' ) '  PDF parameter D = ', d

  if ( .not. burr_check ( a, b, c, d ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call burr_sample ( a, b, c, d, seed, x )

    call burr_pdf ( x, a, b, c, d, pdf )

    call burr_cdf ( x, a, b, c, d, cdf )

    call burr_cdf_inv ( cdf, a, b, c, d, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test027 ( )

!*****************************************************************************80
!
!! TEST027 tests BURR_MEAN, BURR_VARIANCE, BURR_SAMPLE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical burr_check
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST027'
  write ( *, '(a)' ) '  For the Burr PDF:'
  write ( *, '(a)' ) '  BURR_MEAN computes the mean;'
  write ( *, '(a)' ) '  BURR_VARIANCE computes the variance;'
  write ( *, '(a)' ) '  BURR_SAMPLE samples;'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00
  d = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C = ', c
  write ( *, '(a,g14.6)' ) '  PDF parameter D = ', d

  if ( .not. burr_check ( a, b, c, d ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call burr_mean ( a, b, c, d, mean )
  call burr_variance ( a, b, c, d, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call burr_sample ( a, b, c, d, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test0275 ( )

!*****************************************************************************80
!
!! TEST0275 tests CARDIOID_CDF, CARDIOID_CDF_INV, CARDIOID_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 July 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) :: a = 0.0D+00
  real ( kind = 8 ) :: b = 0.25D+00
  logical cardioid_check
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0275'
  write ( *, '(a)' ) '  For the Cardioid PDF:'
  write ( *, '(a)' ) '  CARDIOID_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  CARDIOID_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  CARDIOID_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. cardioid_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call cardioid_sample ( a, b, seed, x )
    call cardioid_pdf ( x, a, b, pdf )
    call cardioid_cdf ( x, a, b, cdf )
    call cardioid_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test0276 ( )

!*****************************************************************************80
!
!! TEST0276 tests CARDIOID_MEAN, CARDIOID_SAMPLE, CARDIOID_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) :: a = 0.0D+00
  real ( kind = 8 ) :: b = 0.25D+00
  logical cardioid_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0276'
  write ( *, '(a)' ) '  For the Cardioid PDF:'
  write ( *, '(a)' ) '  CARDIOID_MEAN computes the mean;'
  write ( *, '(a)' ) '  CARDIOID_SAMPLE samples;'
  write ( *, '(a)' ) '  CARDIOID_VARIANCE computes the variance.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. cardioid_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call cardioid_mean ( a, b, mean )
  call cardioid_variance ( a, b, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call cardioid_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test028 ( )

!*****************************************************************************80
!
!! TEST028 tests CAUCHY_CDF, CAUCHY_CDF_INV, CAUCHY_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical cauchy_check
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST028'
  write ( *, '(a)' ) '  For the Cauchy PDF:'
  write ( *, '(a)' ) '  CAUCHY_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  CAUCHY_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  CAUCHY_PDF evaluates the PDF;'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. cauchy_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call cauchy_sample ( a, b, seed, x )

    call cauchy_pdf ( x, a, b, pdf )

    call cauchy_cdf ( x, a, b, cdf )

    call cauchy_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test029 ( )

!*****************************************************************************80
!
!! TEST029 tests CAUCHY_MEAN, CAUCHY_SAMPLE, CAUCHY_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical cauchy_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST029'
  write ( *, '(a)' ) '  For the Cauchy PDF:'
  write ( *, '(a)' ) '  CAUCHY_MEAN computes the mean;'
  write ( *, '(a)' ) '  CAUCHY_VARIANCE computes the variance;'
  write ( *, '(a)' ) '  CAUCHY_SAMPLE samples.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. cauchy_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call cauchy_mean ( a, b, mean )
  call cauchy_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF mean =        ', variance

  do i = 1, sample_num
    call cauchy_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test030 ( )

!*****************************************************************************80
!
!! TEST030 tests CHI_CDF, CHI_CDF_INV, CHI_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  logical chi_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST030'
  write ( *, '(a)' ) '  For the Chi PDF:'
  write ( *, '(a)' ) '  CHI_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  CHI_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  CHI_PDF evaluates the PDF.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. chi_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call chi_sample ( a, b, c, seed, x )

    call chi_pdf ( x, a, b, c, pdf )

    call chi_cdf ( x, a, b, c, cdf )

    call chi_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test031 ( )

!*****************************************************************************80
!
!! TEST031 tests CHI_MEAN, CHI_SAMPLE, CHI_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical chi_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST031'
  write ( *, '(a)' ) '  For the Chi PDF:'
  write ( *, '(a)' ) '  CHI_MEAN computes the mean;'
  write ( *, '(a)' ) '  CHI_VARIANCE computes the variance;'
  write ( *, '(a)' ) '  CHI_SAMPLE samples.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C = ', c

  if ( .not. chi_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call chi_mean ( a, b, c, mean )
  call chi_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call chi_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test032 ( )

!*****************************************************************************80
!
!! TEST032 tests CHI_SQUARE_CDF, CHI_SQUARE_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) a2
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST032:'
  write ( *, '(a)' ) '  CHI_SQUARE_CDF evaluates the cumulative'
  write ( *, '(a)' ) '    distribution function for the chi-square central'
  write ( *, '(a)' ) '    probability density function.'
  write ( *, '(a)' ) '  CHI_SQUARE_CDF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A    X       Exact F     CHI_SQUARE_CDF(A,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    a2 = real ( a, kind = 8 )

    call chi_square_cdf ( x, a2, fx2 )

    write ( *, '(2x,i4,f8.4,2g14.6)' ) a, x, fx, fx2

  end do

  return
end
subroutine test033 ( )

!*****************************************************************************80
!
!! TEST033 tests CHI_SQUARE_CDF, CHI_SQUARE_CDF_INV, CHI_SQUARE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  logical chi_square_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST033'
  write ( *, '(a)' ) '  For the central chi square PDF:'
  write ( *, '(a)' ) '  CHI_SQUARE_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  CHI_SQUARE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  CHI_SQUARE_PDF evaluates the PDF;'

  a = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. chi_square_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call chi_square_sample ( a, seed, x )

    call chi_square_pdf ( x, a, pdf )

    call chi_square_cdf ( x, a, cdf )

    call chi_square_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test034 ( )

!*****************************************************************************80
!
!! TEST034 tests CHI_SQUARE_MEAN, CHI_SQUARE_SAMPLE, CHI_SQUARE_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  logical chi_square_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST034'
  write ( *, '(a)' ) '  For the central chi square PDF:'
  write ( *, '(a)' ) '  CHI_SQUARE_MEAN computes the mean;'
  write ( *, '(a)' ) '  CHI_SQUARE_SAMPLE samples;'
  write ( *, '(a)' ) '  CHI_SQUARE_VARIANCE computes the variance.'

  a = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. chi_square_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call chi_square_mean ( a, mean )
  call chi_square_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call chi_square_sample ( a, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test035 ( )

!*****************************************************************************80
!
!! TEST035 tests CHI_SQUARE_NONCENTRAL_MEAN, CHI_SQUARE_NONCENTRAL_SAMPLE, CHI_SQUARE_NONCENTRAL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical chi_square_noncentral_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST035'
  write ( *, '(a)' ) '  For the noncentral chi square PDF:'
  write ( *, '(a)' ) '  CHI_SQUARE_NONCENTRAL_SAMPLE samples.'

  a = 3.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. chi_square_noncentral_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call chi_square_noncentral_mean ( a, b, mean )
  call chi_square_noncentral_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)'   ) '  Initial seed =    ', seed

  do i = 1, sample_num
    call chi_square_noncentral_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a,i12)'   ) '  Final seed =      ', seed
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test036 ( )

!*****************************************************************************80
!
!! TEST036 tests CIRCLE_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean(2)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance(2)
  real ( kind = 8 ) x_table(sample_num,2)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xmax(2)
  real ( kind = 8 ) xmin(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST036'
  write ( *, '(a)' ) '  CIRCLE_SAMPLE samples points in a circle.'

  a = 10.0D+00
  b = 4.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X coordinate of center is A = ', a
  write ( *, '(a,g14.6)' ) '  Y coordinate of center is B = ', b
  write ( *, '(a,g14.6)' ) '  Radius is C =                 ', c

  do i = 1, sample_num
    call circle_sample ( a, b, c, seed, x1, x2 )
    x_table(i,1) = x1
    x_table(i,2) = x2
  end do

  do j = 1, 2
    call r8vec_mean ( sample_num, x_table(1,j), mean(j) )
    call r8vec_variance ( sample_num, x_table(1,j), variance(j) )
    call r8vec_max ( sample_num, x_table(1,j), imax, xmax(j) )
    call r8vec_min ( sample_num, x_table(1,j), imin, xmin(j) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'     ) '  Sample size =     ', sample_num
  write ( *, '(a,2g14.6)' ) '  Sample mean =     ', mean(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample variance = ', variance(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample maximum =  ', xmax(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample minimum =  ', xmin(1:2)

  return
end
subroutine test037 ( )

!*****************************************************************************80
!
!! TEST037 tests CIRCULAR_NORMAL_01_*.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean(2)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance(2)
  real ( kind = 8 ) x(2)
  real ( kind = 8 ) x_table(sample_num,2)
  real ( kind = 8 ) xmax(2)
  real ( kind = 8 ) xmin(2)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST037'
  write ( *, '(a)' ) '  For the Circular Normal 01 PDF:'
  write ( *, '(a)' ) '  CIRCULAR_NORMAL_01_MEAN computes the mean;'
  write ( *, '(a)' ) '  CIRCULAR_NORMAL_01_SAMPLE samples;'
  write ( *, '(a)' ) '  CIRCULAR_NORMAL_01_VARIANCE computes variance.'

  call circular_normal_01_mean ( mean )
  call circular_normal_01_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  PDF means =               ', mean(1:2)
  write ( *, '(a,2g14.6)' ) '  PDF variances =           ', variance(1:2)

  do i = 1, sample_num
    call circular_normal_01_sample ( seed, x )
    x_table(i,1) = x(1)
    x_table(i,2) = x(2)
  end do

  do j = 1, 2
    call r8vec_mean ( sample_num, x_table(1,j), mean(j) )
    call r8vec_variance ( sample_num, x_table(1,j), variance(j) )
    call r8vec_max ( sample_num, x_table(1,j), imax, xmax(j) )
    call r8vec_min ( sample_num, x_table(1,j), imin, xmin(j) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'     ) '  Sample size =     ', sample_num
  write ( *, '(a,2g14.6)' ) '  Sample mean =     ', mean(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample variance = ', variance(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample maximum =  ', xmax(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample minimum =  ', xmin(1:2)

  return
end
subroutine test0375 ( )

!*****************************************************************************80
!
!! TEST0375 tests CIRCULAR_NORMAL*.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    03 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a(2)
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean(2)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance(2)
  real ( kind = 8 ) x(2)
  real ( kind = 8 ) x_table(sample_num,2)
  real ( kind = 8 ) xmax(2)
  real ( kind = 8 ) xmin(2)

  a(1) = 1.0D+00
  a(2) = 5.0D+00
  b = 0.75D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0375'
  write ( *, '(a)' ) '  For the Circular Normal PDF:'
  write ( *, '(a)' ) '  CIRCULAR_NORMAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  CIRCULAR_NORMAL_SAMPLE samples;'
  write ( *, '(a)' ) '  CIRCULAR_NORMAL_VARIANCE computes variance.'

  call circular_normal_mean ( a, b, mean )
  call circular_normal_variance ( a, b, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,2g14.6)' ) '  PDF means =               ', mean(1:2)
  write ( *, '(a,2g14.6)' ) '  PDF variances =           ', variance(1:2)

  do i = 1, sample_num
    call circular_normal_sample ( a, b, seed, x )
    x_table(i,1) = x(1)
    x_table(i,2) = x(2)
  end do

  do j = 1, 2
    call r8vec_mean ( sample_num, x_table(1,j), mean(j) )
    call r8vec_variance ( sample_num, x_table(1,j), variance(j) )
    call r8vec_max ( sample_num, x_table(1,j), imax, xmax(j) )
    call r8vec_min ( sample_num, x_table(1,j), imin, xmin(j) )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'     ) '  Sample size =     ', sample_num
  write ( *, '(a,2g14.6)' ) '  Sample mean =     ', mean(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample variance = ', variance(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample maximum =  ', xmax(1:2)
  write ( *, '(a,2g14.6)' ) '  Sample minimum =  ', xmin(1:2)

  return
end
subroutine test038 ( )

!*****************************************************************************80
!
!! TEST038 tests COSINE_CDF, COSINE_CDF_INV, COSINE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  logical cosine_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST038'
  write ( *, '(a)' ) '  For the Cosine PDF:'
  write ( *, '(a)' ) '  COSINE_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  COSINE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  COSINE_PDF evaluates the PDF.'

  a = 2.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. cosine_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call cosine_sample ( a, b, seed, x )

    call cosine_pdf ( x, a, b, pdf )

    call cosine_cdf ( x, a, b, cdf )

    call cosine_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test039 ( )

!*****************************************************************************80
!
!! TEST039 tests COSINE_MEAN, COSINE_SAMPLE, COSINE_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical cosine_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST039'
  write ( *, '(a)' ) '  For the Cosine PDF:'
  write ( *, '(a)' ) '  COSINE_MEAN computes the mean;'
  write ( *, '(a)' ) '  COSINE_SAMPLE samples;'
  write ( *, '(a)' ) '  COSINE_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. cosine_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call cosine_mean ( a, b, mean )
  call cosine_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call cosine_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test0395 ( )

!*****************************************************************************80
!
!! TEST0395 tests COUPON_COMPLETE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) box_num
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) type_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0395'
  write ( *, '(a)' ) '  COUPON_COMPLETE_PDF evaluates the coupon collector''s'
  write ( *, '(a)' ) '  complete collection pdf.'
  write ( *, '(a)' ) ' '

  do type_num = 2, 4

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of coupon types is ', type_num
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '   BOX_NUM      PDF             CDF'
    write ( *, '(a)' ) ' '
    cdf = 0.0D+00
    do box_num = 1, 20
      call coupon_complete_pdf ( type_num, box_num, pdf )
      cdf = cdf + pdf
      write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) box_num, pdf, cdf
    end do

  end do

  return
end
subroutine test040 ( )

!*****************************************************************************80
!
!! TEST040 tests COUPON_SIMULATE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n_trial = 10
  integer ( kind = 4 ), parameter :: max_type = 25

  real ( kind = 8 ) average
  integer ( kind = 4 ) coupon(max_type)
  real ( kind = 8 ) expect
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n_coupon
  integer ( kind = 4 ) n_type
  integer ( kind = 4 ) :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST040'
  write ( *, '(a)' ) '  COUPON_SIMULATE simulates the coupon '
  write ( *, '(a)' ) '  collector''s problem.'
  write ( *, '(a)' ) ' '

  do n_type = 5, max_type, 5

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of coupon types is ', n_type
    expect = real ( n_type, kind = 8 ) * log ( real ( n_type, kind = 8 ) )
    write ( *, '(a,g14.6)' ) '  Expected wait is about ', expect
    write ( *, '(a)' ) ' '

    average = 0.0D+00
    do i = 1, n_trial
      call coupon_simulate ( n_type, seed, coupon, n_coupon )
      write ( *, '(2i5)' ) i, n_coupon
      average = average + real ( n_coupon, kind = 8 )
    end do

    average = average / real ( n_trial, kind = 8 )
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Average wait was ', average

  end do

  return
end
subroutine test041 ( )

!*****************************************************************************80
!
!! TEST041 tests DERANGED_CDF, DERANGED_CDF_INV and DERANGED_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) cdf
  logical deranged_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST041'
  write ( *, '(a)' ) '  For the Deranged PDF:'
  write ( *, '(a)' ) '  DERANGED_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  DERANGED_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  DERANGED_PDF evaluates the PDF;'

  a = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A = ', a

  if ( .not. deranged_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call deranged_sample ( a, seed, x )

    call deranged_pdf ( x, a, pdf )

    call deranged_cdf ( x, a, cdf )

    call deranged_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test042 ( )

!*****************************************************************************80
!
!! TEST042 tests DERANGED_CDF and DERANGED_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  real ( kind = 8 ) cdf
  logical deranged_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  a = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST042'
  write ( *, '(a)' ) '  For the Deranged PDF:'
  write ( *, '(a)' ) '  DERANGED_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  DERANGED_CDF evaluates the CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PDF parameter A = ', a

  if ( .not. deranged_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X    PDF(X)      CDF(X)'
  write ( *, '(a)' ) ' '

  do x = 0, a
    call deranged_pdf ( x, a, pdf )
    call deranged_cdf ( x, a, cdf )
    write ( *, '(2x,i8,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test043 ( )

!*****************************************************************************80
!
!! TEST043 tests DERANGED_MEAN, DERANGED_VARIANCE and DERANGED_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) a
  logical deranged_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST043'
  write ( *, '(a)' ) '  For the Deranged PDF:'
  write ( *, '(a)' ) '  DERANGED_MEAN computes the mean.'
  write ( *, '(a)' ) '  DERANGED_VARIANCE computes the variance.'
  write ( *, '(a)' ) '  DERANGED_SAMPLE samples.'

  a = 7

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A =             ', a

  if ( .not. deranged_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call deranged_mean ( a, mean )
  call deranged_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call deranged_sample ( a, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test044 ( )

!*****************************************************************************80
!
!! TEST044 tests DIGAMMA and PSI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) digamma
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST044:'
  write ( *, '(a)' ) '  DIGAMMA evaluates the DIGAMMA or PSI function.'
  write ( *, '(a)' ) '  PSI_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       DIGAMMA(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    if ( x <= 0.0D+00 ) then
      cycle
    end if

    fx2 = digamma ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test045 ( )

!*****************************************************************************80
!
!! TEST045 tests DIPOLE_CDF, DIPOLE_CDF_INV and DIPOLE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ) atest(test_num)
  real ( kind = 8 ) b
  real ( kind = 8 ) btest(test_num)
  real ( kind = 8 ) cdf
  real ( kind = 8 ) r8_pi
  logical dipole_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) itest
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST045'
  write ( *, '(a)' ) '  For the Dipole PDF:'
  write ( *, '(a)' ) '  DIPOLE_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  DIPOLE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  DIPOLE_PDF evaluates the PDF.'

  atest(1) = 0.0D+00
  btest(1) = 1.0D+00
  atest(2) = r8_pi() / 4.0D+00
  btest(2) = 0.5D+00
  atest(3) = r8_pi() / 2.0D+00
  btest(3) = 0.0D+00

  do itest = 1, test_num

    a = atest(itest)
    b = btest(itest)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
    write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

    if ( .not. dipole_check ( a, b ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Fatal error!'
      write ( *, '(a)' ) '  The parameters are not legal.'
      return
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
    write ( *, '(a)' ) ' '

    do i = 1, 10

      call dipole_sample ( a, b, seed, x )

      call dipole_pdf ( x, a, b, pdf )

      call dipole_cdf ( x, a, b, cdf )

      call dipole_cdf_inv ( cdf, a, b, x2 )

      write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

    end do

  end do

  return
end
subroutine test046 ( )

!*****************************************************************************80
!
!! TEST046 tests DIPOLE_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 10000
  integer ( kind = 4 ), parameter :: test_num = 3

  real ( kind = 8 ) a
  real ( kind = 8 ), dimension ( test_num ) :: a_test = (/ &
    0.0D+00, 0.785398163397448D+00, 1.57079632679490D+00 /)
  real ( kind = 8 ) b
  real ( kind = 8 ), dimension ( test_num ) :: b_test = (/ &
    1.0D+00, 0.5D+00, 0.0D+00 /)
  real ( kind = 8 ) r8_pi
  logical dipole_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST046'
  write ( *, '(a)' ) '  For the Dipole PDF:'
  write ( *, '(a)' ) '  DIPOLE_SAMPLE samples.'

  do test = 1, test_num

    a = a_test(test)
    b = b_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
    write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

    if ( .not. dipole_check ( a, b ) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Fatal error!'
      write ( *, '(a)' ) '  The parameters are not legal.'
      return
    end if

    do i = 1, sample_num
      call dipole_sample ( a, b, seed, x(i) )
    end do

    call r8vec_mean ( sample_num, x, mean )
    call r8vec_variance ( sample_num, x, variance )
    call r8vec_max ( sample_num, x, imax, xmax )
    call r8vec_min ( sample_num, x, imin, xmin )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
    write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
    write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
    write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
    write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  end do

  return
end
subroutine test047 ( )

!*****************************************************************************80
!
!! TEST047 tests DIRICHLET_MEAN, DIRICHLET_SAMPLE and DIRICHLET_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3
  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a(n)
  logical dirichlet_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax(n)
  integer ( kind = 4 ) imin(n)
  real ( kind = 8 ) mean(n)
  real ( kind = 8 ) m2(n,n)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance(n)
  real ( kind = 8 ) x(n,sample_num)
  real ( kind = 8 ) xmax(n)
  real ( kind = 8 ) xmin(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST047'
  write ( *, '(a)' ) '  For the Dirichlet PDF:'
  write ( *, '(a)' ) '  DIRICHLET_SAMPLE samples;'
  write ( *, '(a)' ) '  DIRICHLET_MEAN computes the mean;'
  write ( *, '(a)' ) '  DIRICHLET_VARIANCE computes the variance.'

  a(1:n) = (/ 0.250D+00, 0.500D+00, 1.250D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components N =        ', n

  call r8vec_print ( n, a, '  PDF parameters A:' )
  write ( *, '(a)'    ) '  PDF parameters A(1:N):'

  if ( .not. dirichlet_check ( n, a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call dirichlet_mean ( n, a, mean )

  call dirichlet_variance ( n, a, variance )

  call r8vec_print ( n, mean, '  PDF mean:' )

  call r8vec_print ( n, variance, '  PDF variance:' )

  call dirichlet_moment2 ( n, a, m2 )

  call r8mat_print ( n, n, m2, '  Second moments:' )

  do i = 1, sample_num
    call dirichlet_sample ( n, a, seed, x(1,i) )
  end do

  call r8row_max ( n, sample_num, x, imax, xmax )
  call r8row_min ( n, sample_num, x, imin, xmin )
  call r8row_mean ( n, sample_num, x, mean )
  call r8row_variance ( n, sample_num, x, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Sample size = ', sample_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Observed Mean, Variance, Max, Min:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,4g14.6)' ) i, mean(i), variance(i), xmax(i), xmin(i)
  end do

  return
end
subroutine test048 ( )

!*****************************************************************************80
!
!! TEST048 tests DIRICHLET_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 3

  real ( kind = 8 ) a(n)
  logical dirichlet_check
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST048'
  write ( *, '(a)' ) '  For the Dirichlet PDF:'
  write ( *, '(a)' ) '  DIRICHLET_PDF evaluates the PDF.'

  a(1:3) = (/ 0.250D+00, 0.500D+00, 1.250D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components N =        ', n

  call r8vec_print ( n, a, '  PDF parameters A:' )

  if ( .not. dirichlet_check ( n, a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  x(1:3) = (/ 0.500D+00, 0.125D+00, 0.375D+00 /)

  call r8vec_print ( n, x, '  PDF argument X: ' )

  call dirichlet_pdf ( x, n, a, pdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF value =           ', pdf

  return
end
subroutine test049 ( )

!*****************************************************************************80
!
!! TEST049 tests DIRICHLET_MIX_MEAN and DIRICHLET_MIX_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: comp_num = 2
  integer ( kind = 4 ), parameter :: elem_num = 3
  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a(elem_num,comp_num)
  integer ( kind = 4 ) comp
  real ( kind = 8 ) comp_weight(comp_num)
  logical dirichlet_mix_check
  integer ( kind = 4 ) elem_i
  integer ( kind = 4 ) imax(elem_num)
  integer ( kind = 4 ) imin(elem_num)
  integer ( kind = 4 ) j
  real ( kind = 8 ) mean(elem_num)
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance(elem_num)
  real ( kind = 8 ) x(elem_num,sample_num)
  real ( kind = 8 ) xmax(elem_num)
  real ( kind = 8 ) xmin(elem_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST049'
  write ( *, '(a)' ) '  For the Dirichlet Mixture PDF:'
  write ( *, '(a)' ) '  DIRICHLET_MIX_SAMPLE samples;'
  write ( *, '(a)' ) '  DIRICHLET_MIX_MEAN computes the mean;'

  a(1,1) = 0.250D+00
  a(2,1) = 0.500D+00
  a(3,1) = 1.250D+00

  a(1,2) = 1.500D+00
  a(2,2) = 0.500D+00
  a(3,2) = 2.000D+00

  comp_weight(1) = 1.0D+00
  comp_weight(2) = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of elements ELEM_NUM =   ', elem_num
  write ( *, '(a,i8)' ) '  Number of components COMP_NUM = ', comp_num
  call r8mat_print ( elem_num, comp_num, a, '  PDF parameters A(ELEM,COMP):' )
  call r8vec_print ( comp_num, comp_weight, '  Component weights' )

  if ( .not. dirichlet_mix_check ( comp_num, elem_num, a, comp_weight ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call dirichlet_mix_mean ( comp_num, elem_num, a, comp_weight, mean )

  call r8vec_print ( elem_num, mean, '  PDF means: ' )

  do j = 1, sample_num
    call dirichlet_mix_sample ( comp_num, elem_num, a, &
      comp_weight, seed, comp, x(1,j) )
  end do

  call r8row_max ( elem_num, sample_num, x, imax, xmax )

  call r8row_min ( elem_num, sample_num, x, imin, xmin )

  call r8row_mean ( elem_num, sample_num, x, mean )

  call r8row_variance ( elem_num, sample_num, x, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Sample size = ', sample_num
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Observed Mean, Variance, Max, Min:'
  write ( *, '(a)' ) ' '

  do elem_i = 1, elem_num
    write ( *, '(2x,i8,4g14.6)' ) elem_i, &
      mean(elem_i), variance(elem_i), xmax(elem_i), xmin(elem_i)
  end do

  return
end
subroutine test050 ( )

!*****************************************************************************80
!
!! TEST050 tests DIRICHLET_MIX_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: comp_num = 2
  integer ( kind = 4 ), parameter :: elem_num = 3

  real ( kind = 8 ) a(elem_num,comp_num)
  real ( kind = 8 ) comp_weight(comp_num)
  logical dirichlet_mix_check
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x(elem_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST050'
  write ( *, '(a)' ) '  For the Dirichlet mixture PDF:'
  write ( *, '(a)' ) '  DIRICHLET_MIX_PDF evaluates the PDF.'

  a(1,1) = 0.250D+00
  a(2,1) = 0.500D+00
  a(3,1) = 1.250D+00

  a(1,2) = 1.500D+00
  a(2,2) = 0.500D+00
  a(3,2) = 2.000D+00

  comp_weight(1:2) = (/ 1.0D+00, 2.0D+00 /)

  write ( *, '(a,i8)' ) '  Number of elements ELEM_NUM =   ', elem_num
  write ( *, '(a,i8)' ) '  Number of components COMP_NUM = ', comp_num
  call r8mat_print ( elem_num, comp_num, a, '  PDF parameters A(ELEM,COMP):' )
  call r8vec_print ( comp_num, comp_weight, '  Component weights' )

  if ( .not. dirichlet_mix_check ( comp_num, elem_num, a, comp_weight ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  x(1:3) = (/ 0.500D+00, 0.125D+00, 0.375D+00 /)

  call r8vec_print ( elem_num, x, '  PDF argument X: ' )

  call dirichlet_mix_pdf ( x, comp_num, elem_num, a, comp_weight, &
    pdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF value =           ', pdf

  return
end
subroutine test051 ( )

!*****************************************************************************80
!
!! TEST051 tests BETA_PDF and DIRICHLET_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 2

  real ( kind = 8 ) a
  real ( kind = 8 ) aval
  real ( kind = 8 ) avec(n)
  real ( kind = 8 ) b
  real ( kind = 8 ) bval
  logical dirichlet_check
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(n)

  xval = 0.25D+00
  aval = 2.50D+00
  bval = 3.50D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST051'
  write ( *, '(a)' ) '  BETA_PDF evaluates the Beta PDF.'
  write ( *, '(a)' ) '  DIRICHLET_PDF evaluates the Dirichlet PDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For N = 2, Dirichlet = Beta.'

  xvec(1) = xval
  xvec(2) = 1.0D+00 - xval

  avec(1:2) = (/ aval, bval /)

  if ( .not. dirichlet_check ( n, avec ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of components N =        ', n

  call r8vec_print ( n, avec, '  PDF parameter A: ' )

  call r8vec_print ( n, x, '  PDF argument X: ' )

  call dirichlet_pdf ( xvec, n, avec, pdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Dirichlet PDF value =  ', pdf

  x = xval

  a = aval
  b = bval

  call beta_pdf ( x, a, b, pdf )

  write ( *, '(a,g14.6)' ) '  Beta PDF value =       ', pdf

  return
end
subroutine test052 ( )

!*****************************************************************************80
!
!! TEST052 tests DISCRETE_CDF, DISCRETE_CDF_INV and DISCRETE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 6

  real ( kind = 8 ) b(a)
  real ( kind = 8 ) cdf
  logical discrete_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST052'
  write ( *, '(a)' ) '  For the Discrete PDF:'
  write ( *, '(a)' ) '  DISCRETE_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  DISCRETE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  DISCRETE_PDF evaluates the PDF;'

  b(1:6) = (/ 1.0D+00, 2.0D+00, 6.0D+00, 2.0D+00, 4.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PDF parameter A = ', a
  call r8vec_print ( a, b, '  PDF parameters B = ' )

  if ( .not. discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call discrete_sample ( a, b, seed, x )

    call discrete_pdf ( x, a, b, pdf )

    call discrete_cdf ( x, a, b, cdf )

    call discrete_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test053 ( )

!*****************************************************************************80
!
!! TEST053 tests DISCRETE_MEAN, DISCRETE_SAMPLE and DISCRETE_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 6
  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) b(a)
  logical discrete_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST053'
  write ( *, '(a)' ) '  For the Discrete PDF:'
  write ( *, '(a)' ) '  DISCRETE_MEAN computes the mean;'
  write ( *, '(a)' ) '  DISCRETE_SAMPLE samples;'
  write ( *, '(a)' ) '  DISCRETE_VARIANCE computes the variance.'

  b(1:6) = (/ 1.0D+00, 2.0D+00, 6.0D+00, 2.0D+00, 4.0D+00, 1.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PDF parameter A =             ', a
  call r8vec_print ( a, b, '  PDF parameters B = ' )

  if ( .not. discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call discrete_mean ( a, b, mean )
  call discrete_variance ( a, b, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call discrete_sample ( a, b, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test054 ( )

!*****************************************************************************80
!
!! TEST054 tests EMPIRICAL_DISCRETE_CDF, EMPIRICAL_DISCRETE_CDF_INV, and EMPIRICAL_DISCRETE_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 6

  real ( kind = 8 ), save, dimension ( a ) :: b = (/ &
    1.0D+00, 1.0D+00, 3.0D+00, 2.0D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = 8 ), save, dimension ( a ) :: c = (/ &
    0.0D+00, 1.0D+00, 2.0D+00, 4.5D+00, 6.0D+00, 10.0D+00 /)
  real ( kind = 8 ) cdf
  logical empirical_discrete_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST054'
  write ( *, '(a)' ) '  For the Empirical Discrete PDF:'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A =             ', a
  call r8vec_print ( a, b, '  PDF parameter B:' )
  call r8vec_print ( a, c, '  PDF parameter C:' )

  if ( .not. empirical_discrete_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call empirical_discrete_sample ( a, b, c, seed, x )

    call empirical_discrete_pdf ( x, a, b, c, pdf )

    call empirical_discrete_cdf ( x, a, b, c, cdf )

    call empirical_discrete_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test055 ( )

!*****************************************************************************80
!
!! TEST055 tests EMPIRICAL_DISCRETE_MEAN, EMPIRICAL_DISCRETE_SAMPLE and EMPIRICAL_DISCRETE_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 6
  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ), save, dimension ( a ) :: b = (/ &
    1.0D+00, 1.0D+00, 3.0D+00, 2.0D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = 8 ), save, dimension ( a ) :: c = (/ &
    0.0D+00, 1.0D+00, 2.0D+00, 4.5D+00, 6.0D+00, 10.0D+00 /)
  logical empirical_discrete_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST055'
  write ( *, '(a)' ) '  For the Empirical Discrete PDF:'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_MEAN computes the mean;'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_SAMPLE samples;'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_VARIANCE computes the variance.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'   ) '  PDF parameter A =             ', a
  call r8vec_print ( a, b, '  PDF parameter B:' )
  call r8vec_print ( a, c, '  PDF parameter C:' )

  if ( .not. empirical_discrete_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call empirical_discrete_mean ( a, b, c, mean )
  call empirical_discrete_variance ( a, b, c, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call empirical_discrete_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test056 ( )

!*****************************************************************************80
!
!! TEST056 tests EMPIRICAL_DISCRETE_CDF and EMPIRICAL_DISCRETE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: a = 6

  real ( kind = 8 ), save, dimension ( a ) :: b = (/ &
    1.0D+00, 1.0D+00, 3.0D+00, 2.0D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = 8 ), save, dimension ( a ) :: c = (/ &
    0.0D+00, 1.0D+00, 2.0D+00, 4.5D+00, 6.0D+00, 10.0D+00 /)
  real ( kind = 8 ) cdf
  logical empirical_discrete_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST056'
  write ( *, '(a)' ) '  For the Empirical Discrete PDF.'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  EMPIRICAL_DISCRETE_CDF evaluates the CDF.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PDF parameter A =             ', a
  call r8vec_print ( a, b, '  PDF parameter B:' )
  call r8vec_print ( a, c, '  PDF parameter C:' )

  if ( .not. empirical_discrete_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      PDF(X)      CDF(X)'
  write ( *, '(a)' ) ' '

  do i = -2, 12
    x = real ( i, kind = 8 )
    call empirical_discrete_pdf ( x, a, b, c, pdf )
    call empirical_discrete_cdf ( x, a, b, c, cdf )
    write ( *, '(2x,f8.4,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test0563 ( )

!*****************************************************************************80
!
!! TEST0563 tests ENGLISH_SENTENCE_LENGTH_CDF, ENGLISH_SENTENCE_LENGTH_CDF_INV and ENGLISH_SENTENCE_LENGTH_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0563'
  write ( *, '(a)' ) '  For the English Sentence Length PDF:'
  write ( *, '(a)' ) '  ENGLISH_SENTENCE_LENGTH_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  ENGLISH_SENTENCE_LENGTH_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call english_sentence_length_sample ( seed, x )

    call english_sentence_length_pdf ( x, pdf )

    call english_sentence_length_cdf ( x, cdf )

    call english_sentence_length_cdf_inv ( cdf, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test0564 ( )

!*****************************************************************************80
!
!! TEST0564 tests ENGLISH_SENTENCE_LENGTH_MEAN, ENGLISH_SENTENCE_LENGTH_SAMPLE and ENGLISH_SENTENCE_LENGTH_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    01 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0564'
  write ( *, '(a)' ) '  For the English Sentence Length PDF:'
  write ( *, '(a)' ) '  ENGLISH_SENTENCE_LENGTH_MEAN computes the mean;'
  write ( *, '(a)' ) '  ENGLISH_SENTENCE_LENGTH_SAMPLE samples;'
  write ( *, '(a)' ) '  ENGLISH_SENTENCE_LENGTH_VARIANCE computes the variance.'

  call english_sentence_length_mean ( mean )
  call english_sentence_length_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call english_sentence_length_sample ( seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test0565 ( )

!*****************************************************************************80
!
!! TEST0565 tests ENGLISH_WORD_LENGTH_CDF, ENGLISH_WORD_LENGTH_CDF_INV and ENGLISH_WORD_LENGTH_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0565'
  write ( *, '(a)' ) '  For the English Word Length PDF:'
  write ( *, '(a)' ) '  ENGLISH_WORD_LENGTH_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  ENGLISH_WORD_LENGTH_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  ENGLISH_WORD_LENGTH_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call english_word_length_sample ( seed, x )

    call english_word_length_pdf ( x, pdf )

    call english_word_length_cdf ( x, cdf )

    call english_word_length_cdf_inv ( cdf, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test0566 ( )

!*****************************************************************************80
!
!! TEST0566 tests ENGLISH_WORD_LENGTH_MEAN, ENGLISH_WORD_LENGTH_SAMPLE and ENGLISH_WORD_LENGTH_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0566'
  write ( *, '(a)' ) '  For the English Word Length PDF:'
  write ( *, '(a)' ) '  ENGLISH_WORD_LENGTH_MEAN computes the mean;'
  write ( *, '(a)' ) '  ENGLISH_WORD_LENGTH_SAMPLE samples;'
  write ( *, '(a)' ) '  ENGLISH_WORD_LENGTH_VARIANCE computes the variance.'

  call english_word_length_mean ( mean )
  call english_word_length_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call english_word_length_sample ( seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test057 ( )

!*****************************************************************************80
!
!! TEST057 tests ERLANG_CDF, ERLANG_CDF_INV and ERLANG_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  real ( kind = 8 ) cdf
  logical erlang_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST057'
  write ( *, '(a)' ) '  For the Erlang PDF:'
  write ( *, '(a)' ) '  ERLANG_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  ERLANG_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  ERLANG_PDF evaluates the PDF.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,i8)'    ) '  PDF parameter C = ', c

  if ( .not. erlang_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call erlang_sample ( a, b, c, seed, x )

    call erlang_pdf ( x, a, b, c, pdf )

    call erlang_cdf ( x, a, b, c, cdf )

    call erlang_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test058 ( )

!*****************************************************************************80
!
!! TEST058 tests ERLANG_MEAN, ERLANG_SAMPLE and ERLANG_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) c
  logical erlang_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST058'
  write ( *, '(a)' ) '  For the Erlang PDF:'
  write ( *, '(a)' ) '  ERLANG_MEAN computes the mean;'
  write ( *, '(a)' ) '  ERLANG_SAMPLE samples;'
  write ( *, '(a)' ) '  ERLANG_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =         ', b
  write ( *, '(a,i8)'    ) '  PDF parameter C =         ', c

  if ( .not. erlang_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call erlang_mean ( a, b, c, mean )
  call erlang_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =              ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =          ', variance

  do i = 1, sample_num
    call erlang_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test059 ( )

!*****************************************************************************80
!
!! TEST059 tests ERROR_F and ERROR_F_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error_f
  real ( kind = 8 ) error_f_inverse
  integer ( kind = 4 ) i
  real ( kind = 8 ) r8_normal_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST059'
  write ( *, '(a)' ) '  ERROR_F evaluates the error function erf(x).'
  write ( *, '(a)' ) '  ERROR_F_INVERSE inverts the error function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'X   -> Y = error_F(X) -> Z = error_f_inverse(Y)'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 20

    call normal_01_sample ( seed, x )
    y = error_f ( x )
    z = error_f_inverse ( y )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, z

  end do

  return
end
subroutine test060 ( )

!*****************************************************************************80
!
!! TEST060 tests EXPONENTIAL_01_CDF, EXPONENTIAL_01_CDF_INV, EXPONENTIAL_01_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST060'
  write ( *, '(a)' ) '  For the Exponential 01 PDF:'
  write ( *, '(a)' ) '  EXPONENTIAL_01_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  EXPONENTIAL_01_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  EXPONENTIAL_01_PDF evaluates the PDF.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call exponential_01_sample ( seed, x )

    call exponential_01_pdf ( x, pdf )

    call exponential_01_cdf ( x, cdf )

    call exponential_01_cdf_inv ( cdf, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test061 ( )

!*****************************************************************************80
!
!! TEST061 tests EXPONENTIAL_01_MEAN, EXPONENTIAL_01_SAMPLE, EXPONENTIAL_01_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) mean
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST061'
  write ( *, '(a)' ) '  For the Exponential 01_PDF:'
  write ( *, '(a)' ) '  EXPONENTIAL_01_MEAN computes the mean;'
  write ( *, '(a)' ) '  EXPONENTIAL_01_SAMPLE samples;'
  write ( *, '(a)' ) '  EXPONENTIAL_01_VARIANCE computes the variance.'

  call exponential_01_mean ( mean )
  call exponential_01_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call exponential_01_sample ( seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test062 ( )

!*****************************************************************************80
!
!! TEST062 tests EXPONENTIAL_CDF, EXPONENTIAL_CDF_INV, EXPONENTIAL_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  logical exponential_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST062'
  write ( *, '(a)' ) '  For the Exponential CDF:'
  write ( *, '(a)' ) '  EXPONENTIAL_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  EXPONENTIAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  EXPONENTIAL_PDF evaluates the PDF.'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =         ', b

  if ( .not. exponential_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call exponential_sample ( a, b, seed, x )

    call exponential_pdf ( x, a, b, pdf )

    call exponential_cdf ( x, a, b, cdf )

    call exponential_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test063 ( )

!*****************************************************************************80
!
!! TEST063 tests EXPONENTIAL_MEAN, EXPONENTIAL_SAMPLE, EXPONENTIAL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical exponential_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST063'
  write ( *, '(a)' ) '  For the Exponential PDF:'
  write ( *, '(a)' ) '  EXPONENTIAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  EXPONENTIAL_SAMPLE samples;'
  write ( *, '(a)' ) '  EXPONENTIAL_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =       ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. exponential_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call exponential_mean ( a, b, mean )
  call exponential_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call exponential_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test064 ( )

!*****************************************************************************80
!
!! TEST064 tests EXTREME_VALUES_CDF, EXTREME_VALUES_CDF_INV, EXTREME_VALUES_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  logical extreme_values_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST064'
  write ( *, '(a)' ) '  For the Extreme Values CDF:'
  write ( *, '(a)' ) '  EXTREME_VALUES_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  EXTREME_VALUES_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  EXTREME_VALUES_PDF evaluates the PDF;'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. extreme_values_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call extreme_values_sample ( a, b, seed, x )

    call extreme_values_pdf ( x, a, b, pdf )

    call extreme_values_cdf ( x, a, b, cdf )

    call extreme_values_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test065 ( )

!*****************************************************************************80
!
!! TEST065 tests EXTREME_VALUES_MEAN, *_SAMPLE, *_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical extreme_values_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST065'
  write ( *, '(a)' ) '  For the Extreme Values PDF:'
  write ( *, '(a)' ) '  EXTREME_VALUES_MEAN computes the mean;'
  write ( *, '(a)' ) '  EXTREME_VALUES_SAMPLE samples;'
  write ( *, '(a)' ) '  EXTREME_VALUES_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. extreme_values_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call extreme_values_mean ( a, b, mean )
  call extreme_values_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  do i = 1, sample_num
    call extreme_values_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test066 ( )

!*****************************************************************************80
!
!! TEST066 tests F_CDF and F_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST066:'
  write ( *, '(a)' ) '  F_CDF evaluates the F central CDF.'
  write ( *, '(a)' ) '  F_CDF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     A       B    X       Exact F       F_CDF(A,B,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call f_cdf ( x, a, b, fx2 )

    write ( *, '(2x,i8,2x,i8,2x,f8.4,2g14.6)' ) a, b, x, fx, fx2

  end do

  return
end
subroutine test067 ( )

!*****************************************************************************80
!
!! TEST067 tests F_CDF, F_PDF and F_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  logical f_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST067'
  write ( *, '(a)' ) '  For the central F PDF:'
  write ( *, '(a)' ) '  F_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  F_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  F_SAMPLE samples the PDF.'

  m = 1
  n = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter M = ', m
  write ( *, '(a,i8)'    ) '  PDF parameter N = ', n

  if ( .not. f_check ( m, n ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call f_sample ( m, n, seed, x )

    call f_pdf ( x, m, n, pdf )

    call f_cdf ( x, m, n, cdf )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test068 ( )

!*****************************************************************************80
!
!! TEST068 tests F_MEAN, F_SAMPLE, F_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  logical f_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) m
  real ( kind = 8 ) mean
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST068'
  write ( *, '(a)' ) '  For the central F PDF:'
  write ( *, '(a)' ) '  F_MEAN computes the mean;'
  write ( *, '(a)' ) '  F_SAMPLE samples;'
  write ( *, '(a)' ) '  F_VARIANCE computes the varianc.'

  m = 8
  n = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter M =             ', m
  write ( *, '(a,i8)'    ) '  PDF parameter N =             ', n

  if ( .not. f_check ( m, n ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call f_mean ( m, n, mean )
  call f_variance ( m, n, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =          ', variance

  do i = 1, sample_num
    call f_sample ( m, n, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test069 ( )

!*****************************************************************************80
!
!! TEST069 tests FACTORIAL_LOG and GAMMA_LOG_INT;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) f
  real ( kind = 8 ) factorial_log
  real ( kind = 8 ) g
  real ( kind = 8 ) gamma_log_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_1 = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST069'
  write ( *, '(a)' ) &
    '  FACTORIAL_LOG evaluates the log of the factorial function;'
  write ( *, '(a)' ) '  GAMMA_LOG_INT evaluates the log for integer argument.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I    GAMMA_LOG_INT(I+1)   FACTORIAL_LOG(I)'
  write ( *, '(a)' ) ' '

  do i = 1, 20
    g = gamma_log_int ( i+i4_1 )
    f = factorial_log ( i )
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, g, f
  end do

  return
end
subroutine test070 ( )

!*****************************************************************************80
!
!! TEST070 tests FACTORIAL_STIRLING and I4_FACTORIAL;
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

  integer ( kind = 4 ) i4_factorial
  real ( kind = 8 ) factorial_stirling
  integer ( kind = 4 ) i
  real ( kind = 8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST070'
  write ( *, '(a)' ) '  FACTORIAL_STIRLING computes Stirling''s'
  write ( *, '(a)' ) '    approximate factorial function;'
  write ( *, '(a)' ) '  I4_FACTORIAL evaluates the factorial function;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N        Stirling       N!'
  write ( *, '(a)' ) ' '

  do i = 0, 20
    value = factorial_stirling ( i )
    write ( *, '(2x,i8,2x,g14.6,2x,i20)' ) i, value, i4_factorial ( i )
  end do

  return
end
subroutine test07025 ( )

!*****************************************************************************80
!
!! TEST07025 tests FERMI_DIRAC_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    25 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 10000
  integer ( kind = 4 ), parameter :: test_num = 7

  integer ( kind = 4 ) i
  real ( kind = 8 ) mean
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  real ( kind = 8 ) u
  real ( kind = 8 ), dimension ( test_num ) :: u_test = (/ &
   1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00, 16.0D+00, &
  32.0D+00, 1.0D+00  /)
  real ( kind = 8 ) v
  real ( kind = 8 ), dimension ( test_num ) :: v_test = (/ &
   1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00, &
   1.0D+00, 0.25D+00  /)
  real ( kind = 8 ) variance
  real ( kind = 8 ) z(sample_num)
  real ( kind = 8 ) z_max
  real ( kind = 8 ) z_min

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07025'
  write ( *, '(a)' ) '  Test FERMI_DIRAC_SAMPLE:'

  do test = 1, test_num

    u = u_test(test)
    v = v_test(test)
    seed = 123456789

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  U =          ', u
    write ( *, '(a,g14.6)' ) '  V =          ', v
    write ( *, '(a,i8)'    ) '  SAMPLE_NUM = ', sample_num
    write ( *, '(a,i12)'   ) '  SEED =       ', seed

    do i = 1, sample_num
      call fermi_dirac_sample ( u, v, seed, z(i) )
    end do

    z_max = maxval ( z(1:sample_num) )
    z_min = minval ( z(1:sample_num) )

    call r8vec_mean ( sample_num, z, mean )
    call r8vec_variance ( sample_num, z, variance )

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Minimum value =   ', z_min
    write ( *, '(a,g14.6)' ) '  Maximum value =   ', z_max
    write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
    write ( *, '(a,g14.6)' ) '  Sample variance = ', variance

  end do

  return
end
subroutine test0705 ( )

!*****************************************************************************80
!
!! TEST0705 tests FISHER_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: n = 10
  integer ( kind = 4 ), parameter :: test_num = 3

  integer ( kind = 4 ) j
  real ( kind = 8 ) kappa
  real ( kind = 8 ) mu(3)
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) seed
  integer ( kind = 4 ) test
  real ( kind = 8 ) x(3,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0705'
  write ( *, '(a)' ) '  For the Fisher PDF:'
  write ( *, '(a)' ) '  FISHER_SAMPLE samples the PDF.'
  write ( *, '(a)' ) '  FISHER_PDF evaluates the PDF.'

  do test = 1, test_num

    if ( test == 1 ) then
      kappa = 0.0D+00
      mu = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
    else if ( test == 2 ) then
      kappa = 0.5D+00
      mu = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
    else if ( test == 3 ) then
      kappa = 10.0D+00
      mu = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  PDF parameters:'
    write ( *, '(a,g14.6)' ) '    Concentration parameter KAPPA =      ', kappa
    write ( *, '(a,3f8.4)' ) '    Direction MU(1:3) = ', mu(1:3)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '      X                         PDF'
    write ( *, '(a)' ) ' '

    seed = 123456789

    call fisher_sample ( kappa, mu, n, seed, x )

    do j = 1, n

      call fisher_pdf ( x(1:3,j), kappa, mu, pdf )

      write ( *, '(2x,3f8.4,2x,g14.6)' ) x(1:3,j), pdf

    end do

  end do

  return
end
subroutine test071 ( )

!*****************************************************************************80
!
!! TEST071 tests FISK_CDF, FISK_CDF_INV and FISK_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  logical fisk_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST071'
  write ( *, '(a)' ) '  For the Fisk PDF:'
  write ( *, '(a)' ) '  FISK_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  FISK_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  FISK_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. fisk_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call fisk_sample ( a, b, c, seed, x )

    call fisk_pdf ( x, a, b, c, pdf )

    call fisk_cdf ( x, a, b, c, cdf )

    call fisk_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test072 ( )

!*****************************************************************************80
!
!! TEST072 tests FISK_MEAN, FISK_SAMPLE and FISK_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical fisk_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST072'
  write ( *, '(a)' ) '  For the Fisk PDF:'
  write ( *, '(a)' ) '  FISK_MEAN computes the mean;'
  write ( *, '(a)' ) '  FISK_SAMPLE samples;'
  write ( *, '(a)' ) '  FISK_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. fisk_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call fisk_mean ( a, b, c, mean )
  call fisk_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call fisk_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test073 ( )

!*****************************************************************************80
!
!! TEST073 tests FOLDED_NORMAL_CDF, FOLDED_NORMAL_CDF_INV, FOLDED_NORMAL_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  logical folded_normal_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST073'
  write ( *, '(a)' ) '  For the Folded Normal PDF:'
  write ( *, '(a)' ) '  FOLDED_NORMAL_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  FOLDED_NORMAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  FOLDED_NORMAL_PDF evaluates the PDF.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =         ', b

  if ( .not. folded_normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call folded_normal_sample ( a, b, seed, x )

    call folded_normal_pdf ( x, a, b, pdf )

    call folded_normal_cdf ( x, a, b, cdf )

    call folded_normal_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test074 ( )

!*****************************************************************************80
!
!! TEST074 tests FOLDED_NORMAL_MEAN, FOLDED_NORMAL_SAMPLE, FOLDED_NORMAL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical folded_normal_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST074'
  write ( *, '(a)' ) '  For the Folded Normal PDF:'
  write ( *, '(a)' ) '  FOLDED_NORMAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  FOLDED_NORMAL_SAMPLE samples;'
  write ( *, '(a)' ) '  FOLDED_NORMAL_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. folded_normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call folded_normal_mean ( a, b, mean )
  call folded_normal_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call folded_normal_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test0744 ( )

!*****************************************************************************80
!
!! TEST0744 tests FRECHET_CDF, FRECHET_CDF_INV and FRECHET_PDF.
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
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) alpha
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0744'
  write ( *, '(a)' ) '  For the Frechet PDF:'
  write ( *, '(a)' ) '  FRECHET_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  FRECHET_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  FRECHET_PDF evaluates the PDF;'

  alpha = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter ALPHA =         ', alpha

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call frechet_sample ( alpha, seed, x )

    call frechet_pdf ( x, alpha, pdf )

    call frechet_cdf ( x, alpha, cdf )

    call frechet_cdf_inv ( cdf, alpha, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test0745 ( )

!*****************************************************************************80
!
!! TEST0745 tests FRECHET_MEAN, FRECHET_SAMPLE and FRECHET_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 September 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) alpha
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0745'
  write ( *, '(a)' ) '  For the Frechet PDF:'
  write ( *, '(a)' ) '  FRECHET_MEAN computes the mean;'
  write ( *, '(a)' ) '  FRECHET_SAMPLE samples;'
  write ( *, '(a)' ) '  FRECHET_VARIANCE computes the variance.'

  alpha = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter ALPHA =         ', alpha

  call frechet_mean ( alpha, mean )
  call frechet_variance ( alpha, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call frechet_sample ( alpha, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test075 ( )

!*****************************************************************************80
!
!! TEST075 tests GAMMA, GAMMA_LOG, GAMMA_LOG_INT, I_FACTORIAL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) i4_factorial
  real ( kind = 8 ) g1
  real ( kind = 8 ) g2
  real ( kind = 8 ) g3
  real ( kind = 8 ) g4
  real ( kind = 8 ) gamma
  real ( kind = 8 ) gamma_log
  real ( kind = 8 ) gamma_log_int
  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_1 = 1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST075'
  write ( *, '(a)' ) '  GAMMA evaluates the Gamma function;'
  write ( *, '(a)' ) '  GAMMA_LOG evaluates the log of the Gamma function;'
  write ( *, '(a)' ) '  GAMMA_LOG_INT evaluates the log for integer argument;'
  write ( *, '(a)' ) '  I_FACTORIAL evaluates the factorial function.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  X, GAMMA(X), Exp(GAMMA_LOG(X)), Exp(GAMMA_LOG_INT(X)) ' // &
    'I_FACTORIAL(X+1)'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    x = real ( i, kind = 8 )
    g1 = gamma ( x )
    g2 = exp ( gamma_log ( x ) )
    g3 = exp ( gamma_log_int ( i ) )
    g4 = i4_factorial ( i - i4_1 )
    write ( *, '(2x,5g14.6)' ) x, g1, g2, g3, g4
  end do

  return
end
subroutine test076 ( )

!*****************************************************************************80
!
!! TEST076 tests GAMMA_INC and GAMMA_INC_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  real ( kind = 8 ) gamma_inc
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST076:'
  write ( *, '(a)' ) '  GAMMA_INC evaluates the normalized incomplete Gamma'
  write ( *, '(a)' ) '    function P(A,X).'
  write ( *, '(a)' ) '  GAMMA_INC_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A     X       Exact F      GAMMA_INC(A,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = gamma_inc ( a, x )

    write ( *, '(2x,2f8.4,2g14.6)' ) a, x, fx, fx2

  end do

  return
end
subroutine test077 ( )

!*****************************************************************************80
!
!! TEST077 tests GAMMA_CDF, GAMMA_PDF, GAMMA_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  logical              gamma_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST077'
  write ( *, '(a)' ) '  For the Gamma PDF:'
  write ( *, '(a)' ) '  GAMMA_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  GAMMA_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  GAMMA_SAMPLE samples the PDF.'

  a = 1.0D+00
  b = 1.5D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. gamma_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X  PDF   CDF'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call gamma_sample ( a, b, c, seed, x )

    call gamma_cdf ( x, a, b, c, cdf )

    call gamma_pdf ( x, a, b, c, pdf )

    write ( *, '(2x,3g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test078 ( )

!*****************************************************************************80
!
!! TEST078 tests GAMMA_MEAN, GAMMA_SAMPLE and GAMMA_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical              gamma_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST078'
  write ( *, '(a)' ) '  For the Gamma PDF:'
  write ( *, '(a)' ) '  GAMMA_MEAN computes the mean;'
  write ( *, '(a)' ) '  GAMMA_SAMPLE samples;'
  write ( *, '(a)' ) '  GAMMA_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 3.0D+00
  c = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. gamma_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call gamma_mean ( a, b, c, mean )
  call gamma_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call gamma_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test079 ( )

!*****************************************************************************80
!
!! TEST079 tests GENLOGISTIC_CDF, GENLOGISTIC_CDF_INV, GENLOGISTIC_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  logical              genlogistic_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST079'
  write ( *, '(a)' ) '  For the Generalized Logistic PDF:'
  write ( *, '(a)' ) '  GENLOGISTIC_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  GENLOGISTIC_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  GENLOGISTIC_CDF_INV inverts the CDF.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. genlogistic_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call genlogistic_sample ( a, b, c, seed, x )

    call genlogistic_pdf ( x, a, b, c, pdf )

    call genlogistic_cdf ( x, a, b, c, cdf )

    call genlogistic_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test080 ( )

!*****************************************************************************80
!
!! TEST080 tests GENLOGISTIC_MEAN, GENLOGISTIC_SAMPLE, GENLOGISTIC_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical              genlogistic_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST080'
  write ( *, '(a)' ) '  For the Generalized Logistic PDF:'
  write ( *, '(a)' ) '  GENLOGISTIC_MEAN computes the mean;'
  write ( *, '(a)' ) '  GENLOGISTIC_SAMPLE samples;'
  write ( *, '(a)' ) '  GENLOGISTIC_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. genlogistic_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call genlogistic_mean ( a, b, c, mean )
  call genlogistic_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call genlogistic_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test081 ( )

!*****************************************************************************80
!
!! TEST081 tests GEOMETRIC_CDF, GEOMETRIC_CDF_INV, GEOMETRIC_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  logical              geometric_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST081'
  write ( *, '(a)' ) '  For the Geometric PDF:'
  write ( *, '(a)' ) '  GEOMETRIC_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  GEOMETRIC_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  GEOMETRIC_PDF evaluates the PDF;'

  a = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =      ', a

  if ( .not. geometric_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call geometric_sample ( a, seed, x )

    call geometric_pdf ( x, a, pdf )

    call geometric_cdf ( x, a, cdf )

    call geometric_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test082 ( )

!*****************************************************************************80
!
!! TEST082 tests GEOMETRIC_MEAN, GEOMETRIC_SAMPLE, GEOMETRIC_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  logical              geometric_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST082'
  write ( *, '(a)' ) '  For the Geometric PDF:'
  write ( *, '(a)' ) '  GEOMETRIC_MEAN computes the mean;'
  write ( *, '(a)' ) '  GEOMETRIC_SAMPLE samples;'
  write ( *, '(a)' ) '  GEOMETRIC_VARIANCE computes the variance.'

  a = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. geometric_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call geometric_mean ( a, mean )
  call geometric_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call geometric_sample ( a, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test083 ( )

!*****************************************************************************80
!
!! TEST083 tests GEOMETRIC_CDF and GEOMETRIC_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  logical              geometric_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST083'
  write ( *, '(a)' ) '  For the Geometric PDF:'
  write ( *, '(a)' ) '  GEOMETRIC_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  GEOMETRIC_CDF evaluates the CDF.'

  a = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a

  if ( .not. geometric_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      PDF(X)      CDF(X)'
  write ( *, '(a)' ) ' '

  do x = 0, 10
    call geometric_pdf ( x, a, pdf )
    call geometric_cdf ( x, a, cdf )
    write ( *, '(2x,i8,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test084 ( )

!*****************************************************************************80
!
!! TEST084 tests GOMPERTZ_CDF, GOMPERTZ_CDF_INV and GOMPERTZ_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  logical              gompertz_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST084'
  write ( *, '(a)' ) '  For the Gompertz PDF:'
  write ( *, '(a)' ) '  GOMPERTZ_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  GOMPERTZ_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  GOMPERTZ_PDF evaluates the PDF;'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =       ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =       ', b

  if ( .not. gompertz_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call gompertz_sample ( a, b, seed, x )

    call gompertz_pdf ( x, a, b, pdf )

    call gompertz_cdf ( x, a, b, cdf )

    call gompertz_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test085 ( )

!*****************************************************************************80
!
!! TEST085 tests GOMPERTZ_SAMPLE;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical gompertz_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST085'
  write ( *, '(a)' ) '  For the Gompertz PDF:'
  write ( *, '(a)' ) '  GOMPERTZ_SAMPLE samples;'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =       ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =       ', b

  if ( .not. gompertz_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  do i = 1, sample_num
    call gompertz_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test086 ( )

!*****************************************************************************80
!
!! TEST086 tests GUMBEL_CDF, GUMBEL_CDF_INV, GUMBEL_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST086'
  write ( *, '(a)' ) '  For the Gumbel PDF:'
  write ( *, '(a)' ) '  GUMBEL_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  GUMBEL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  GUMBEL_PDF evaluates the PDF.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call gumbel_sample ( seed, x )

    call gumbel_pdf ( x, pdf )

    call gumbel_cdf ( x, cdf )

    call gumbel_cdf_inv ( cdf, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test087 ( )

!*****************************************************************************80
!
!! TEST087 tests GUMBEL_MEAN, GUMBEL_SAMPLE, GUMBEL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST087'
  write ( *, '(a)' ) '  For the Gumbel PDF:'
  write ( *, '(a)' ) '  GUMBEL_MEAN computes the mean;'
  write ( *, '(a)' ) '  GUMBEL_SAMPLE samples;'
  write ( *, '(a)' ) '  GUMBEL_VARIANCE computes the variance.'

  call gumbel_mean ( mean )

  call gumbel_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean      =               ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call gumbel_sample ( seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test088 ( )

!*****************************************************************************80
!
!! TEST088 tests HALF_NORMAL_CDF, HALF_NORMAL_CDF_INV, HALF_NORMAL_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  logical half_normal_check
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST088'
  write ( *, '(a)' ) '  For the Half Normal PDF:'
  write ( *, '(a)' ) '  HALF_NORMAL_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  HALF_NORMAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  HALF_NORMAL_PDF evaluates the PDF.'

  a = 0.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =         ', b

  if ( .not. half_normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call half_normal_sample ( a, b, seed, x )

    call half_normal_pdf ( x, a, b, pdf )

    call half_normal_cdf ( x, a, b, cdf )

    call half_normal_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test089 ( )

!*****************************************************************************80
!
!! TEST089 tests HALF_NORMAL_MEAN, HALF_NORMAL_SAMPLE, HALF_NORMAL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  logical half_normal_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST089'
  write ( *, '(a)' ) '  For the Half Normal PDF:'
  write ( *, '(a)' ) '  HALF_NORMAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  HALF_NORMAL_SAMPLE samples;'
  write ( *, '(a)' ) '  HALF_NORMAL_VARIANCE computes the variance.'

  a = 0.0D+00
  b = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. half_normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call half_normal_mean ( a, b, mean )
  call half_normal_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call half_normal_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test090 ( )

!*****************************************************************************80
!
!! TEST090 tests HYPERGEOMETRIC_CDF and HYPERGEOMETRIC_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) cdf
  logical hypergeometric_check
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST090'
  write ( *, '(a)' ) '  For the Hypergeometric PDF:'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_PDF evaluates the PDF.'

  x = 7

  n = 100
  m = 70
  l = 1000

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Total number of balls =         ', l
  write ( *, '(a,i8)'    ) '  Number of white balls =         ', m
  write ( *, '(a,i8)'    ) '  Number of balls taken =         ', n

  if ( .not. hypergeometric_check ( n, m, l ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call hypergeometric_pdf ( x, n, m, l, pdf )

  call hypergeometric_cdf ( x, n, m, l, cdf )

  write ( *, '(a,i8)'    ) '  PDF argument X =                ', x
  write ( *, '(a,g14.6)' ) '  PDF value =                   = ', pdf
  write ( *, '(a,g14.6)' ) '  CDF value =                   = ', cdf

  return
end
subroutine test091 ( )

!*****************************************************************************80
!
!! TEST091 tests HYPERGEOMETRIC_MEAN, HYPERGEOMETRIC_SAMPLE, HYPERGEOMETRIC_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  logical hypergeometric_check
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ) mean
  integer ( kind = 4 ) n
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST091'
  write ( *, '(a)' ) '  For the Hypergeometric PDF:'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_MEAN computes the mean;'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_SAMPLE samples;'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_VARIANCE computes the variance.'

  n = 100
  m = 70
  l = 1000

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter N =             ', n
  write ( *, '(a,i8)'    ) '  PDF parameter M =             ', m
  write ( *, '(a,i8)'    ) '  PDF parameter L =             ', l

  if ( .not. hypergeometric_check ( n, m, l ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call hypergeometric_mean ( n, m, l, mean )
  call hypergeometric_variance ( n, m, l, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'THIS CALL IS TAKING FOREVER!'
  return

  do i = 1, sample_num
    call hypergeometric_sample ( n, m, l, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test092 ( )

!*****************************************************************************80
!
!! TEST092 tests R8_CEILING.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) r8_ceiling
  integer ( kind = 4 ) ival
  real ( kind = 8 ) rval

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST092'
  write ( *, '(a)' ) '  R8_CEILING rounds an R8 up.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X           R8_CEILING(X)'
  write ( *, '(a)' ) ' '

  do i = -6, 6
    rval = real ( i, kind = 8 ) / 5.0D+00
    ival = r8_ceiling ( rval )
    write ( *, '(2x,g14.6,i8)' ) rval, ival
  end do

  return
end
subroutine test093 ( )

!*****************************************************************************80
!
!! TEST093 tests INVERSE_GAUSSIAN_CDF and INVERSE_GAUSSIAN_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical inverse_gaussian_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST093'
  write ( *, '(a)' ) '  For the Inverse Gaussian PDF:'
  write ( *, '(a)' ) '  INVERSE_GAUSSIAN_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  INVERSE_GAUSSIAN_PDF evaluates the PDF.'

  a = 5.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. inverse_gaussian_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call inverse_gaussian_sample ( a, b, seed, x )

    call inverse_gaussian_pdf ( x, a, b, pdf )

    call inverse_gaussian_cdf ( x, a, b, cdf )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test094 ( )

!*****************************************************************************80
!
!! TEST094 tests INVERSE_GAUSSIAN_MEAN, INVERSE_GAUSSIAN_SAMPLE, INVERSE_GAUSSIAN_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical inverse_gaussian_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST094'
  write ( *, '(a)' ) '  For the Inverse Gaussian PDF:'
  write ( *, '(a)' ) '  INVERSE_GAUSSIAN_MEAN computes the mean;'
  write ( *, '(a)' ) '  INVERSE_GAUSSIAN_SAMPLE samples;'
  write ( *, '(a)' ) '  INVERSE_GAUSSIAN_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. inverse_gaussian_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call inverse_gaussian_mean ( a, b, mean )
  call inverse_gaussian_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call inverse_gaussian_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test095 ( )

!*****************************************************************************80
!
!! TEST095 tests LAPLACE_CDF, LAPLACE_CDF_INV and LAPLACE_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical laplace_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST095'
  write ( *, '(a)' ) '  For the Laplace PDF:'
  write ( *, '(a)' ) '  LAPLACE_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LAPLACE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  LAPLACE_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. laplace_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call laplace_sample ( a, b, seed, x )

    call laplace_pdf ( x, a, b, pdf )

    call laplace_cdf ( x, a, b, cdf )

    call laplace_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test096 ( )

!*****************************************************************************80
!
!! TEST096 tests LAPLACE_MEAN, LAPLACE_SAMPLE, LAPLACE_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical laplace_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST096'
  write ( *, '(a)' ) '  For the Laplace PDF:'
  write ( *, '(a)' ) '  LAPLACE_MEAN computes the mean;'
  write ( *, '(a)' ) '  LAPLACE_SAMPLE samples;'
  write ( *, '(a)' ) '  LAPLACE_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. laplace_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call laplace_mean ( a, b, mean )
  call laplace_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call laplace_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test0965 ( )

!*****************************************************************************80
!
!! TEST0965 tests LEVY_CDF, LEVY_CDF_INV and LEVY_PDF.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) test
  integer ( kind = 4 ), parameter :: test_num = 10
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST0965'
  write ( *, '(a)' ) '  For the Levy PDF:'
  write ( *, '(a)' ) '  LEVY_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LEVY_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  LEVY_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X              PDF           CDF          X2'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    call levy_sample ( a, b, seed, x )

    call levy_pdf ( x, a, b, pdf )

    call levy_cdf ( x, a, b, cdf )

    call levy_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,g12.6,2x,g12.6,2x,g12.6,2x,g12.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test097 ( )

!*****************************************************************************80
!
!! TEST097 tests LOGISTIC_CDF, LOGISTIC_CDF_INV, LOGISTIC_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical logistic_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST097'
  write ( *, '(a)' ) '  For the Logistic PDF:'
  write ( *, '(a)' ) '  LOGISTIC_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LOGISTIC_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  LOGISTIC_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. logistic_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call logistic_sample ( a, b, seed, x )

    call logistic_pdf ( x, a, b, pdf )

    call logistic_cdf ( x, a, b, cdf )

    call logistic_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test098 ( )

!*****************************************************************************80
!
!! TEST098 tests LOGISTIC_MEAN, LOGISTIC_SAMPLE, LOGISTIC_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical logistic_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST098'
  write ( *, '(a)' ) '  For the Logistic PDF:'
  write ( *, '(a)' ) '  LOGISTIC_MEAN computes the mean;'
  write ( *, '(a)' ) '  LOGISTIC_SAMPLE samples;'
  write ( *, '(a)' ) '  LOGISTIC_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. logistic_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call logistic_mean ( a, b, mean )
  call logistic_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call logistic_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test099 ( )

!*****************************************************************************80
!
!! TEST099 tests LOG_NORMAL_CDF, LOG_NORMAL_CDF_INV, LOG_NORMAL_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical log_normal_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST099'
  write ( *, '(a)' ) '  For the Lognormal PDF:'
  write ( *, '(a)' ) '  LOG_NORMAL_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LOG_NORMAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  LOG_NORMAL_PDF evaluates the PDF;'

  a = 10.0D+00
  b = 2.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. log_normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call log_normal_sample ( a, b, seed, x )

    call log_normal_pdf ( x, a, b, pdf )

    call log_normal_cdf ( x, a, b, cdf )

    call log_normal_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test100 ( )

!*****************************************************************************80
!
!! TEST100 tests LOG_NORMAL_MEAN, LOG_NORMAL_SAMPLE, LOG_NORMAL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical log_normal_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST100'
  write ( *, '(a)' ) '  For the Lognormal PDF:'
  write ( *, '(a)' ) '  LOG_NORMAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  LOG_NORMAL_SAMPLE samples;'
  write ( *, '(a)' ) '  LOG_NORMAL_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. log_normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST100 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call log_normal_mean ( a, b, mean )
  call log_normal_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call log_normal_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test101 ( )

!*****************************************************************************80
!
!! TEST101 tests LOG_SERIES_CDF, LOG_SERIES_CDF_INV, LOG_SERIES_PDF.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical log_series_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST101'
  write ( *, '(a)' ) '  For the Logseries PDF,'
  write ( *, '(a)' ) '  LOG_SERIES_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LOG_SERIES_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  LOG_SERIES_PDF evaluates the PDF;'

  a = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =  ', a

  if ( .not. log_series_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST101 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call log_series_sample ( a, seed, x )

    call log_series_pdf ( x, a, pdf )

    call log_series_cdf ( x, a, cdf )

    call log_series_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test102 ( )

!*****************************************************************************80
!
!! TEST102 tests LOG_SERIES_CDF and LOG_SERIES_PDF.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  logical log_series_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST102'
  write ( *, '(a)' ) '  For the Logseries PDF:'
  write ( *, '(a)' ) '  LOG_SERIES_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LOG_SERIES_PDF evaluates the PDF.'

  x = 2

  a = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A  = ', a

  if ( .not. log_series_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST102 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X   PDF(X)       CDF(X)'
  write ( *, '(a)' ) ' '

  do x = 1, 10
    call log_series_pdf ( x, a, pdf )
    call log_series_cdf ( x, a, cdf )
    write ( *, '(2x,i8,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test103 ( )

!*****************************************************************************80
!
!! TEST103 tests LOG_SERIES_MEAN, LOG_SERIES_SAMPLE and LOG_SERIES_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical log_series_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST103'
  write ( *, '(a)' ) '  For the Logseries PDF:'
  write ( *, '(a)' ) '  LOG_SERIES_MEAN computes the mean;'
  write ( *, '(a)' ) '  LOG_SERIES_VARIANCE computes the variance;'
  write ( *, '(a)' ) '  LOG_SERIES_SAMPLE samples.'

  a = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. log_series_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST103 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call log_series_mean ( a, mean )
  call log_series_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call log_series_sample ( a, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test104 ( )

!*****************************************************************************80
!
!! TEST104 tests LOG_UNIFORM_CDF, LOG_UNIFORM_INV, LOG_UNIFORM_PDF;
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical log_uniform_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST104'
  write ( *, '(a)' ) '  For the Log Uniform PDF:'
  write ( *, '(a)' ) '  LOG_UNIFORM_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LOG_UNIFORM_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  LOG_UNIFORM_PDF evaluates the PDF;'

  a = 2.0D+00
  b = 20.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. log_uniform_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST104 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call log_uniform_sample ( a, b, seed, x )

    call log_uniform_pdf ( x, a, b, pdf )

    call log_uniform_cdf ( x, a, b, cdf )

    call log_uniform_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test105 ( )

!*****************************************************************************80
!
!! TEST105 tests LOG_UNIFORM_MEAN and LOG_UNIFORM_SAMPLE;
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical log_uniform_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST105'
  write ( *, '(a)' ) '  For the Log Uniform PDF:'
  write ( *, '(a)' ) '  LOG_UNIFORM_MEAN computes the mean;'
  write ( *, '(a)' ) '  LOG_UNIFORM_SAMPLE samples;'

  a = 2.0D+00
  b = 20.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. log_uniform_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST105 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call log_uniform_mean ( a, b, mean )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean

  do i = 1, sample_num
    call log_uniform_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test106 ( )

!*****************************************************************************80
!
!! TEST106 tests LORENTZ_CDF, LORENTZ_CDF_INV and LORENTZ_PDF.
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

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST106'
  write ( *, '(a)' ) '  For the Lorentz PDF:'
  write ( *, '(a)' ) '  LORENTZ_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  LORENTZ_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  LORENTZ_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call lorentz_sample ( seed, x )

    call lorentz_pdf ( x, pdf )

    call lorentz_cdf ( x, cdf )

    call lorentz_cdf_inv ( cdf, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test107 ( )

!*****************************************************************************80
!
!! TEST107 tests LORENTZ_MEAN, LORENTZ_SAMPLE and LORENTZ_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST107'
  write ( *, '(a)' ) '  For the Lorentz PDF:'
  write ( *, '(a)' ) '  LORENTZ_MEAN computes the mean;'
  write ( *, '(a)' ) '  LORENTZ_VARIANCE computes the variance;'
  write ( *, '(a)' ) '  LORENTZ_SAMPLE samples.'

  call lorentz_mean ( mean )
  call lorentz_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call lorentz_sample ( seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test108 ( )

!*****************************************************************************80
!
!! TEST108 tests MAXWELL_CDF, MAXWELL_CDF_INV and MAXWELL_PDF.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical maxwell_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST108'
  write ( *, '(a)' ) '  For the Maxwell CDF:'
  write ( *, '(a)' ) '  MAXWELL_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  MAXWELL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  MAXWELL_PDF evaluates the PDF.'

  a = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a

  if ( .not. maxwell_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST108 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call maxwell_sample ( a, seed, x )

    call maxwell_pdf ( x, a, pdf )

    call maxwell_cdf ( x, a, cdf )

    call maxwell_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test109 ( )

!*****************************************************************************80
!
!! TEST109 tests MAXWELL_MEAN, MAXWELL_SAMPLE, MAXWELL_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical maxwell_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST109'
  write ( *, '(a)' ) '  For the Maxwell PDF:'
  write ( *, '(a)' ) '  MAXWELL_MEAN computes the mean;'
  write ( *, '(a)' ) '  MAXWELL_VARIANCE computes the variance;'
  write ( *, '(a)' ) '  MAXWELL_SAMPLE samples.'

  a = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. maxwell_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST109 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call maxwell_mean ( a, mean )
  call maxwell_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', variance

  do i = 1, sample_num
    call maxwell_sample ( a, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test110 ( )

!*****************************************************************************80
!
!! TEST110 tests MULTINOMIAL_COEF1, MULTINOMIAL_COEF2.
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

  integer ( kind = 4 ), parameter :: maxfactor = 5

  integer ( kind = 4 ) factor(maxfactor)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncomb1
  integer ( kind = 4 ) ncomb2
  integer ( kind = 4 ) nfactor

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST110'
  write ( *, '(a)' ) '  MULTINOMIAL_COEF1 computes multinomial'
  write ( *, '(a)' ) '    coefficients using the Gamma function;'
  write ( *, '(a)' ) '  MULTINOMIAL_COEF2 computes multinomial'
  write ( *, '(a)' ) '    coefficients directly.'

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

    write ( *, '(i4,i4,3x,i5,i5)' ) factor(1), factor(2), ncomb1, ncomb2

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

      write ( *, '(i4,i4,i4,3x,i5,i5)' ) factor(1), factor(2), factor(3), &
        ncomb1, ncomb2

    end do

  end do

  return
end
subroutine test111 ( )

!*****************************************************************************80
!
!! TEST111 tests MULTINOMIAL_MEAN, MULTINOMIAL_SAMPLE, MULTINOMIAL_VARIANCE;
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

  integer ( kind = 4 ), parameter :: b = 3
  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax(b)
  integer ( kind = 4 ) imin(b)
  real ( kind = 8 ) mean(b)
  logical multinomial_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance(b)
  integer ( kind = 4 ) x(b,sample_num)
  integer ( kind = 4 ) xmax(b)
  integer ( kind = 4 ) xmin(b)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST111'
  write ( *, '(a)' ) '  For the Multinomial PDF:'
  write ( *, '(a)' ) '  MULTINOMIAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  MULTINOMIAL_SAMPLE samples;'
  write ( *, '(a)' ) '  MULTINOMIAL_VARIANCE computes the variance;'

  a = 5

  c(1:3) = (/ 0.125D+00, 0.500D+00, 0.375D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PDF parameter A =             ', a
  write ( *, '(a,i8)' ) '  PDF parameter B =             ', b
  call r8vec_print ( b, c, '  PDF parameter C = ' )

  if ( .not. multinomial_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST111 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call multinomial_mean ( a, b, c, mean )
  call multinomial_variance ( a, b, c, variance )

  call r8vec_print ( b, mean, '  PDF means: ' )
  call r8vec_print ( b, variance, '  PDF variances:' )

  do i = 1, sample_num
    call multinomial_sample ( a, b, c, seed, x(1,i) )
  end do

  call i4row_max ( b, sample_num, x, imax, xmax )
  call i4row_min ( b, sample_num, x, imin, xmin )
  call i4row_mean ( b, sample_num, x, mean )
  call i4row_variance ( b, sample_num, x, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Sample size = ', sample_num
  write ( *, '(a)' ) '  Component Mean, Variance, Min, Max:'
  do i = 1, b
    write ( *, '(2x,i8,2g14.6,2i8)' ) i, mean(i), variance(i), xmin(i), xmax(i)
  end do

  return
end
subroutine test112 ( )

!*****************************************************************************80
!
!! TEST112 tests MULTINOMIAL_PDF.
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

  integer ( kind = 4 ), parameter :: b = 3

  integer ( kind = 4 ) a
  real ( kind = 8 ) c(b)
  integer ( kind = 4 ) i
  logical multinomial_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x(b)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST112'
  write ( *, '(a)' ) '  For the Multinomial PDF:'
  write ( *, '(a)' ) '  MULTINOMIAL_PDF evaluates the PDF.'
  a = 5

  c(1:3) = (/ 0.10D+00, 0.50D+00, 0.40D+00 /)

  write ( *, '(a,i8)' ) '  PDF parameter A = ', a
  write ( *, '(a,i8)' ) '  PDF parameter B = ', b
  call r8vec_print ( b, c, '  PDF parameter C:' )

  if ( .not. multinomial_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST112 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  x(1:3) = (/ 0, 2, 3 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  PDF argument X:'
  write ( *, '(a)' ) ' '
  do i = 1, b
    write ( *, '(2x,i8)' ) x(i)
  end do

  call multinomial_pdf ( x, a, b, c, pdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF value =     ', pdf

  return
end
subroutine test113 ( )

!*****************************************************************************80
!
!! TEST113 tests NAKAGAMI_CDF, NAKAGAMI_PDF.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  logical nakagami_check
  real ( kind = 8 ) pdf
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST113'
  write ( *, '(a)' ) '  For the Nakagami PDF:'
  write ( *, '(a)' ) '  NAKAGAMI_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  NAKAGAMI_PDF evaluates the PDF;'

  x = 1.25D+00

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. nakagami_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST113 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call nakagami_pdf ( x, a, b, c, pdf )

  call nakagami_cdf ( x, a, b, c, cdf )

  write ( *, '(a,g14.6)' ) '  PDF argument X =              ', x
  write ( *, '(a,g14.6)' ) '  PDF value =                   ', pdf
  write ( *, '(a,g14.6)' ) '  CDF value =                   ', cdf

  return
end
subroutine test114 ( )

!*****************************************************************************80
!
!! TEST114 tests NAKAGAMI_MEAN, NAKAGAMI_VARIANCE.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  real ( kind = 8 ) mean
  logical nakagami_check
  real ( kind = 8 ) variance

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST114'
  write ( *, '(a)' ) '  For the Nakagami PDF:'
  write ( *, '(a)' ) '  NAKAGAMI_MEAN computes the mean;'
  write ( *, '(a)' ) '  NAKAGAMI_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. nakagami_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST114 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call nakagami_mean ( a, b, c, mean )
  call nakagami_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  return
end
subroutine test1145 ( )

!*****************************************************************************80
!
!! TEST1145 tests NEGATIVE_BINOMIAL_CDF, NEGATIVE_BINOMIAL_CDF_INV, NEGATIVE_BINOMIAL_PDF.
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

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical negative_binomial_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1145'
  write ( *, '(a)' ) '  For the Negative Binomial PDF:'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_PDF evaluates the PDF.'

  a = 2
  b = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. negative_binomial_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST1145 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call negative_binomial_sample ( a, b, seed, x )

    call negative_binomial_pdf ( x, a, b, pdf )

    call negative_binomial_cdf ( x, a, b, cdf )

    call negative_binomial_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test1146 ( )

!*****************************************************************************80
!
!! TEST1146 tests NEGATIVE_BINOMIAL_MEAN, NEGATIVE_BINOMIAL_SAMPLE, NEGATIVE_BINOMIAL_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical negative_binomial_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1146'
  write ( *, '(a)' ) '  For the Negative Binomial PDF:'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_SAMPLE samples;'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_VARIANCE computes the variance.'

  a = 2
  b = 0.75D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. negative_binomial_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST1146 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call negative_binomial_mean ( a, b, mean )
  call negative_binomial_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call negative_binomial_sample ( a, b, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test115 ( )

!*****************************************************************************80
!
!! TEST115 tests NORMAL_01_CDF, NORMAL_01_CDF_INV, NORMAL_01_PDF;
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

  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST115'
  write ( *, '(a)' ) '  For the Normal 01 PDF:'
  write ( *, '(a)' ) '  NORMAL_01_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  NORMAL_01_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  NORMAL_01_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call normal_01_sample ( seed, x )

    call normal_01_pdf ( x, pdf )

    call normal_01_cdf ( x, cdf )

    call normal_01_cdf_inv ( cdf, x2 )

    write ( *, '(2x,g24.16,2x,g14.6,2x,g14.6,2x,g24.16)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test116 ( )

!*****************************************************************************80
!
!! TEST116 tests NORMAL_01_MEAN, NORMAL_01_SAMPLE, NORMAL_01_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST116'
  write ( *, '(a)' ) '  For the Normal 01 PDF:'
  write ( *, '(a)' ) '  NORMAL_01_MEAN computes the mean;'
  write ( *, '(a)' ) '  NORMAL_01_SAMPLE samples the PDF;'
  write ( *, '(a)' ) '  NORMAL_01_VARIANCE returns the variance.'

  call normal_01_mean ( mean )
  call normal_01_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call normal_01_sample ( seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test117 ( )

!*****************************************************************************80
!
!! TEST117 tests NORMAL_CDF, NORMAL_CDF_INV, NORMAL_PDF;
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical normal_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST117'
  write ( *, '(a)' ) '  For the Normal PDF:'
  write ( *, '(a)' ) '  NORMAL_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  NORMAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  NORMAL_PDF evaluates the PDF;'

  a = 100.0D+00
  b = 15.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST117 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call normal_sample ( a, b, seed, x )

    call normal_pdf ( x, a, b, pdf )

    call normal_cdf ( x, a, b, cdf )

    call normal_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test118 ( )

!*****************************************************************************80
!
!! TEST118 tests NORMAL_MEAN, NORMAL_SAMPLE, NORMAL_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  logical normal_check
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST118'
  write ( *, '(a)' ) '  For the Normal PDF:'
  write ( *, '(a)' ) '  NORMAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  NORMAL_SAMPLE samples;'
  write ( *, '(a)' ) '  NORMAL_VARIANCE returns the variance.'

  a = 100.0D+00
  b = 15.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. normal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST118 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call normal_mean ( a, b, mean )
  call normal_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =        ', variance

  do i = 1, sample_num
    call normal_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test119 ( )

!*****************************************************************************80
!
!! TEST119 tests PARETO_CDF, PARETO_CDF_INV, PARETO_PDF.
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
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  logical pareto_check
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST119'
  write ( *, '(a)' ) '  For the Pareto PDF:'
  write ( *, '(a)' ) '  PARETO_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  PARETO_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  PARETO_PDF evaluates the PDF;'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. pareto_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST119 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call pareto_sample ( a, b, seed, x )

    call pareto_pdf ( x, a, b, pdf )

    call pareto_cdf ( x, a, b, cdf )

    call pareto_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test120 ( )

!*****************************************************************************80
!
!! TEST120 tests PARETO_MEAN, PARETO_SAMPLE, PARETO_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical pareto_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST120'
  write ( *, '(a)' ) '  For the Pareto PDF:'
  write ( *, '(a)' ) '  PARETO_MEAN computes the mean;'
  write ( *, '(a)' ) '  PARETO_SAMPLE samples;'
  write ( *, '(a)' ) '  PARETO_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. pareto_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST120 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call pareto_mean ( a, b, mean )
  call pareto_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call pareto_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test123

!*****************************************************************************80
!
!! TEST123 tests PEARSON_05_PDF.
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
  real ( kind = 8 ) c
  real ( kind = 8 ) pdf
  logical pearson_05_check
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST123'
  write ( *, '(a)' ) '  For the Pearson 05 PDF:'
  write ( *, '(a)' ) '  PEARSON_05_PDF evaluates the PDF.'

  x = 5.0D+00

  a = 1.0D+00
  b = 2.0D+00
  c = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C = ', c

  if ( .not. pearson_05_check ( a, b, c ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST123 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call pearson_05_pdf ( x, a, b, c, pdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF argument X =  ', x
  write ( *, '(a,g14.6)' ) '  PDF value =       ', pdf

  return
end
subroutine test124

!*****************************************************************************80
!
!! TEST124 tests PLANCK_PDF, PLANCK_SAMPLE.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  logical planck_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST124'
  write ( *, '(a)' ) '  For the Planck PDF:'
  write ( *, '(a)' ) '  PLANCK_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  PLANCK_SAMPLE samples the PDF.'

  a = 2.0D+00;
  b = 3.0D+00;

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. planck_check ( a, b ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST124 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call planck_sample ( a, b, seed, x )

    call planck_pdf ( x, a, b, pdf )

    write ( *, '(2x,2g14.6)' ) x, pdf

  end do

  return
end
subroutine test125

!*****************************************************************************80
!
!! TEST125 tests PLANCK_MEAN, PLANCK_SAMPLE, PLANCK_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical planck_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST125'
  write ( *, '(a)' ) '  For the Planck PDF:'
  write ( *, '(a)' ) '  PLANCK_MEAN computes the mean.'
  write ( *, '(a)' ) '  PLANCK_SAMPLE samples.'
  write ( *, '(a)' ) '  PLANCK_VARIANCE computes the variance.'
  write ( *, '(a)' ) ' '

  a = 2.0D+00;
  b = 3.0D+00;

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. planck_check ( a, b ) ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST125 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call planck_mean ( a, b, mean )
  call planck_variance ( a, b, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =     ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance = ', variance

  do i = 1, sample_num
    call planck_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test126

!*****************************************************************************80
!
!! TEST126 tests POISSON_CDF, POISSON_CDF_VALUES.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST126:'
  write ( *, '(a)' ) '  POISSON_CDF evaluates the cumulative distribution'
  write ( *, '(a)' ) '    function for the discrete Poisson probability'
  write ( *, '(a)' ) '    density function.'
  write ( *, '(a)' ) '  POISSON_CDF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A is the expected mean number of successes per unit time;'
  write ( *, '(a)' ) '  X is the number of successes;'
  write ( *, '(a)' ) '  POISSON_CDF is the probability of having up to X'
  write ( *, '(a)' ) '  successes in unit time.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   A          X   Exact F     POISSON_CDF(A,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call poisson_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call poisson_cdf ( x, a, fx2 )

    write ( *, '(2x,f8.4,i8,2g14.6)' ) a, x, fx, fx2

  end do

  return
end
subroutine test127

!*****************************************************************************80
!
!! TEST127 tests POISSON_CDF, POISSON_CDF_INV, POISSON_PDF.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  logical poisson_check
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST127'
  write ( *, '(a)' ) '  For the Poisson PDF:'
  write ( *, '(a)' ) '  POISSON_CDF evaluates the CDF,'
  write ( *, '(a)' ) '  POISSON_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  POISSON_PDF evaluates the PDF.'

  a = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =  ', a

  if ( .not. poisson_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST127 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call poisson_sample ( a, seed, x )

    call poisson_pdf ( x, a, pdf )

    call poisson_cdf ( x, a, cdf )

    call poisson_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test128

!*****************************************************************************80
!
!! TEST128 tests POISSON_MEAN, POISSON_SAMPLE, POISSON_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical poisson_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST128'
  write ( *, '(a)' ) '  For the Poisson PDF:'
  write ( *, '(a)' ) '  POISSON_MEAN computes the mean;'
  write ( *, '(a)' ) '  POISSON_SAMPLE samples;'
  write ( *, '(a)' ) '  POISSON_VARIANCE computes the variance.'

  a = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. poisson_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST128 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call poisson_mean ( a, mean )
  call poisson_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call poisson_sample ( a, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test129

!*****************************************************************************80
!
!! TEST129 tests POWER_CDF, POWER_CDF_INV, POWER_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  logical power_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST129'
  write ( *, '(a)' ) '  For the Power PDF:'
  write ( *, '(a)' ) '  POWER_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  POWER_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  POWER_PDF evaluates the PDF;'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =       ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =       ', b

  if ( .not. power_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST129 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call power_sample ( a, b, seed, x )

    call power_pdf ( x, a, b, pdf )

    call power_cdf ( x, a, b, cdf )

    call power_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test130

!*****************************************************************************80
!
!! TEST130 tests POWER_MEAN, POWER_SAMPLE, POWER_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical power_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST130'
  write ( *, '(a)' ) '  For the Power PDF:'
  write ( *, '(a)' ) '  POWER_MEAN computes the mean;'
  write ( *, '(a)' ) '  POWER_SAMPLE samples;'
  write ( *, '(a)' ) '  POWER_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. power_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST130 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call power_mean ( a, b, mean )
  call power_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call power_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test1304

!*****************************************************************************80
!
!! TEST1304 tests QUASIGEOMETRIC_CDF, *_CDF_INV, *_PDF;
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  logical              quasigeometric_check
  integer ( kind = 4 ) :: seed = 123456789
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1304'
  write ( *, '(a)' ) '  For the Quasigeometric PDF:'
  write ( *, '(a)' ) '  QUASIGEOMETRIC_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  QUASIGEOMETRIC_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  QUASIGEOMETRIC_PDF evaluates the PDF;'

  a = 0.4825D+00
  b = 0.5893D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =      ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =      ', b

  if ( .not. quasigeometric_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call quasigeometric_sample ( a, b, seed, x )

    call quasigeometric_pdf ( x, a, b, pdf )

    call quasigeometric_cdf ( x, a, b, cdf )

    call quasigeometric_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test1306

!*****************************************************************************80
!
!! TEST1306 tests QUASIGEOMETRIC_MEAN, *_SAMPLE, *_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical              quasigeometric_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1306'
  write ( *, '(a)' ) '  For the Quasigeometric PDF:'
  write ( *, '(a)' ) '  QUASIGEOMETRIC_MEAN computes the mean;'
  write ( *, '(a)' ) '  QUASIGEOMETRIC_SAMPLE samples;'
  write ( *, '(a)' ) '  QUASIGEOMETRIC_VARIANCE computes the variance.'

  a = 0.4825D+00
  b = 0.5893D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =      ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =      ', b

  if ( .not. quasigeometric_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call quasigeometric_mean ( a, b, mean )
  call quasigeometric_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call quasigeometric_sample ( a, b, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test131

!*****************************************************************************80
!
!! TEST131 tests RAYLEIGH_CDF, RAYLEIGH_CDF_INV, RAYLEIGH_PDF.
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

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  logical rayleigh_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST131'
  write ( *, '(a)' ) '  For the Rayleigh PDF:'
  write ( *, '(a)' ) '  RAYLEIGH_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  RAYLEIGH_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  RAYLEIGH_PDF evaluates the PDF;'

  a = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. rayleigh_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST131 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call rayleigh_sample ( a, seed, x )

    call rayleigh_pdf ( x, a, pdf )

    call rayleigh_cdf ( x, a, cdf )

    call rayleigh_cdf_inv ( cdf, a, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test132

!*****************************************************************************80
!
!! TEST132 tests RAYLEIGH_MEAN, RAYLEIGH_SAMPLE, RAYLEIGH_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical rayleigh_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST132'
  write ( *, '(a)' ) '  For the Rayleigh PDF:'
  write ( *, '(a)' ) '  RAYLEIGH_MEAN computes the mean;'
  write ( *, '(a)' ) '  RAYLEIGH_SAMPLE samples;'
  write ( *, '(a)' ) '  RAYLEIGH_VARIANCE computes the variance.'

  a = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. rayleigh_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST132 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call rayleigh_mean ( a, mean )
  call rayleigh_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call rayleigh_sample ( a, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test133

!*****************************************************************************80
!
!! TEST133 tests RECIPROCAL_CDF, RECIPROCAL_CDF_INV, RECIPROCAL_CDF;
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

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  logical reciprocal_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST133'
  write ( *, '(a)' ) '  For the Reciprocal PDF:'
  write ( *, '(a)' ) '  RECIPROCAL_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  RECIPROCAL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  RECIPROCAL_PDF evaluates the PDF.'

  a = 1.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =         ', b

  if ( .not. reciprocal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST133 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call reciprocal_sample ( a, b, seed, x )

    call reciprocal_pdf ( x, a, b, pdf )

    call reciprocal_cdf ( x, a, b, cdf )

    call reciprocal_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test134

!*****************************************************************************80
!
!! TEST134 tests RECIPROCAL_MEAN, RECIPROCAL_SAMPLE, RECIPROCAL_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical reciprocal_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST134'
  write ( *, '(a)' ) '  For the Reciprocal PDF:'
  write ( *, '(a)' ) '  RECIPROCAL_MEAN computes the mean;'
  write ( *, '(a)' ) '  RECIPROCAL_SAMPLE samples;'
  write ( *, '(a)' ) '  RECIPROCAL_VARIANCE computes the variance.'

  a = 1.0D+00
  b = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B = ', b

  if ( .not. reciprocal_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST134 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call reciprocal_mean ( a, b, mean )
  call reciprocal_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call reciprocal_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test1341

!*****************************************************************************80
!
!! TEST1341 checks RIBESL against BESSEL_IX_VALUES.
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

  integer ( kind = 4 ), parameter :: nb_max = 10

  real ( kind = 8 ) alpha
  real ( kind = 8 ) alpha_frac
  real ( kind = 8 ) b(nb_max)
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1341:'
  write ( *, '(a)' ) '  RIBESL computes values of Bessel functions'
  write ( *, '(a)' ) '  of NONINTEGER order.'
  write ( *, '(a)' ) '  BESSEL_IX_VALUES returns selected values of the'
  write ( *, '(a)' ) '  Bessel function In for NONINTEGER order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      ALPHA         X             FX' // &
    '                        FX2'
  write ( *, '(a)' ) '                                  (table)' // &
    '                   (RIBESL)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_ix_values ( n_data, alpha, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    ize = 1

    nb = int ( alpha ) + 1

    if ( nb_max < nb ) then
      write ( *, * ) '  [Skipping calculation, NB_MAX too small.]'
      cycle
    end if

    alpha_frac = alpha - real ( int ( alpha ), kind = 8 )

    call ribesl ( x, alpha_frac, nb, ize, b, ncalc )

    fx2 = b(nb)

    write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16)' ) alpha, x, fx, fx2

  end do

  return
end
subroutine test1342

!*****************************************************************************80
!
!! TEST1342 tests RUNS_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    24 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) pdf
  real ( kind = 8 ) pdf_total
  integer ( kind = 4 ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1342'
  write ( *, '(a)' ) '  For the RUNS PDF:'
  write ( *, '(a)' ) '  RUNS_PDF evaluates the PDF;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  M is the number of symbols of one kind,'
  write ( *, '(a)' ) '  N is the number of symbols of the other kind,'
  write ( *, '(a)' ) '  R is the number of runs (sequences of one symbol)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M         N         R      PDF'
  write ( *, '(a)' ) ' '

  m = 6

  do n = 0, 8

    write ( *, '(a)' ) ' '
    pdf_total = 0.0D+00

    do r = 1, 2 * min ( m, n ) + 2

      call runs_pdf ( m, n, r, pdf )
      write ( *, '(2x,i8,2x,i8,2x,i8,2x,g14.6)' ) m, n, r, pdf
      pdf_total = pdf_total + pdf

    end do

    write ( *, '(2x,i8,2x,8x,2x,8x,2x,g14.6)' ) m, pdf_total

  end do

  return
end
subroutine test1344

!*****************************************************************************80
!
!! TEST1344 tests RUNS_MEAN, RUNS_VARIANCE;
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) m
  real ( kind = 8 ) mean
  integer ( kind = 4 ) n
  integer ( kind = 4 ) r(sample_num)
  integer ( kind = 4 ) rmax
  integer ( kind = 4 ) rmin
  integer ( kind = 4 ) seed
  real ( kind = 8 ) variance

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1344'
  write ( *, '(a)' ) '  For the RUNS PDF:'
  write ( *, '(a)' ) '  RUNS_MEAN computes the mean;'
  write ( *, '(a)' ) '  RUNS_VARIANCE computes the variance'

  m = 10
  n = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter M = ', m
  write ( *, '(a,g14.6)' ) '  PDF parameter N = ', n

  call runs_mean ( m, n, mean )
  call runs_variance ( m, n, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =        ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =    ', variance

  seed = 123456789

  do i = 1, sample_num
    call runs_sample ( m, n, seed, r(i) )
  end do

  call i4vec_mean ( sample_num, r, mean )
  call i4vec_variance ( sample_num, r, variance )
  call i4vec_max ( sample_num, r, imax, rmax )
  call i4vec_min ( sample_num, r, imin, rmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', rmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', rmin


  return
end
subroutine test135

!*****************************************************************************80
!
!! TEST135 tests SECH_CDF, SECH_CDF_INV, SECH_PDF.
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
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  logical sech_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST135'
  write ( *, '(a)' ) '  For the Hyperbolic Secant PDF:'
  write ( *, '(a)' ) '  SECH_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  SECH_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  SECH_PDF evaluates the PDF.'

  a = 3.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =         ', b

  if ( .not. sech_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST135 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call sech_sample ( a, b, seed, x )

    call sech_pdf ( x, a, b, pdf )

    call sech_cdf ( x, a, b, cdf )

    call sech_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test136

!*****************************************************************************80
!
!! TEST136 tests SECH_MEAN, SECH_SAMPLE, SECH_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  logical sech_check
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST136'
  write ( *, '(a)' ) '  For the Hyperbolic Secant PDF:'
  write ( *, '(a)' ) '  SECH_MEAN computes the mean;'
  write ( *, '(a)' ) '  SECH_SAMPLE samples;'
  write ( *, '(a)' ) '  SECH_VARIANCE computes the variance.'

  a = 3.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. sech_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST136 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call sech_mean ( a, b, mean )
  call sech_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call sech_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test137

!*****************************************************************************80
!
!! TEST137 tests SEMICIRCULAR_CDF, SEMICIRCULAR_CDF_INV, SEMICIRCULAR_PDF.
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
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical semicircular_check
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST137'
  write ( *, '(a)' ) '  For the Semicircular PDF:'
  write ( *, '(a)' ) '  SEMICIRCULAR_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  SEMICIRCULAR_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  SEMICIRCULAR_PDF evaluates the PDF.'

  a = 3.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =         ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =         ', b

  if ( .not. semicircular_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST137 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call semicircular_sample ( a, b, seed, x )

    call semicircular_pdf ( x, a, b, pdf )

    call semicircular_cdf ( x, a, b, cdf )

    call semicircular_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test138

!*****************************************************************************80
!
!! TEST138 tests SEMICIRCULAR_MEAN, SEMICIRCULAR_SAMPLE, SEMICIRCULAR_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  logical semicircular_check
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST138'
  write ( *, '(a)' ) '  For the Semicircular PDF:'
  write ( *, '(a)' ) '  SEMICIRCULAR_MEAN computes the mean;'
  write ( *, '(a)' ) '  SEMICIRCULAR_SAMPLE samples;'
  write ( *, '(a)' ) '  SEMICIRCULAR_VARIANCE computes the variance.'

  a = 3.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. semicircular_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST138 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call semicircular_mean ( a, b, mean )
  call semicircular_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call semicircular_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test139

!*****************************************************************************80
!
!! TEST139 tests STUDENT_CDF, STUDENT_CDF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 November 2005
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
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST139:'
  write ( *, '(a)' ) '  STUDENT_CDF evaluates the cumulative'
  write ( *, '(a)' ) '    distribution function for the Student''s central T'
  write ( *, '(a)' ) '    probability density function.'
  write ( *, '(a)' ) '  STUDENT_CDF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   A       B       C       X       Exact F ' // &
    '    STUDENT_CDF(A,B,C,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_cdf_values ( n_data, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    a = 0.0D+00
    b = 1.0D+00

    call student_cdf ( x, a, b, c, fx2 )

    write ( *, '(4f8.4,2g14.6)' ) a, b, c, x, fx, fx2

  end do

  return
end
subroutine test140

!*****************************************************************************80
!
!! TEST140 tests STUDENT_CDF, STUDENT_PDF and STUDENT_SAMPLE.
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
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) seed
  logical student_check
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST140'
  write ( *, '(a)' ) '  For the central Student PDF:'
  write ( *, '(a)' ) '  STUDENT_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  STUDENT_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  STUDENT_SAMPLE samples the PDF.'

  a = 0.5D+00
  b = 2.0D+00
  c = 6.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =   ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =   ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =   ', c

  if ( .not. student_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST140 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF'
  write ( *, '(a)' ) ' '

  seed = 123456789

  do i = 1, 10

    call student_sample ( a, b, c, seed, x )

    call student_pdf ( x, a, b, c, pdf )

    call student_cdf ( x, a, b, c, cdf )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test141

!*****************************************************************************80
!
!! TEST141 tests STUDENT_MEAN, STUDENT_SAMPLE, STUDENT_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  logical student_check
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST141'
  write ( *, '(a)' ) '  For the central Student PDF:'
  write ( *, '(a)' ) '  STUDENT_MEAN computes the mean;'
  write ( *, '(a)' ) '  STUDENT_SAMPLE samples;'
  write ( *, '(a)' ) '  STUDENT_VARIANCE computes the variance.'

  a = 0.5D+00
  b = 2.0D+00
  c = 6.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. student_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST141 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call student_mean ( a, b, c, mean )
  call student_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call student_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test142

!*****************************************************************************80
!
!! TEST142 tests STUDENT_NONCENTRAL_CDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) idf
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST142'
  write ( *, '(a)' ) '  For the Noncentral Student PDF:'
  write ( *, '(a)' ) '  STUDENT_NONCENTRAL_CDF evaluates the CDF;'

  x = 0.50D+00

  idf = 10
  b = 1.0D+00

  call student_noncentral_cdf ( x, idf, b, cdf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF argument X =              ', x
  write ( *, '(a,i8)' ) '  PDF parameter IDF =           ', idf
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  CDF value =                   ', cdf

  return
end
subroutine test1425

!*****************************************************************************80
!
!! TEST1425 tests TFN and OWEN_VALUES.
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
  real ( kind = 8 ) h
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) t
  real ( kind = 8 ) t2
  real ( kind = 8 ) tfn

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1425'
  write ( *, '(a)' ) '  TFN evaluates Owen''s T function;'
  write ( *, '(a)' ) '  OWEN_VALUES stores some exact values.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      H             A           T(H,A)      Exact'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t )

    if ( n_data <= 0 ) then
      exit
    end if

    t2 = tfn ( h, a )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) h, a, t2, t

  end do

  return
end
subroutine test143

!*****************************************************************************80
!
!! TEST143 tests TRIANGLE_CDF, TRIANGLE_CDF_INV and TRIANGLE_PDF;
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
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical triangle_check
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST143'
  write ( *, '(a)' ) '  For the Triangle PDF:'
  write ( *, '(a)' ) '  TRIANGLE_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  TRIANGLE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  TRIANGLE_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 3.0D+00
  c = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =      ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =      ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =      ', c

  if ( .not. triangle_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST143 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call triangle_sample ( a, b, c, seed, x )

    call triangle_pdf ( x, a, b, c, pdf )

    call triangle_cdf ( x, a, b, c, cdf )

    call triangle_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test144

!*****************************************************************************80
!
!! TEST144 tests TRIANGLE_MEAN, TRIANGLE_SAMPLE and TRIANGLE_VARIANCE;
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  logical triangle_check
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST144'
  write ( *, '(a)' ) '  For the Triangle PDF:'
  write ( *, '(a)' ) '  TRIANGLE_MEAN returns the mean;'
  write ( *, '(a)' ) '  TRIANGLE_SAMPLE samples;'
  write ( *, '(a)' ) '  TRIANGLE_VARIANCE returns the variance;'

  a = 1.0D+00
  b = 3.0D+00
  c = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. triangle_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST144 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call triangle_mean ( a, b, c, mean )
  call triangle_variance ( a, b, c, variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter MEAN =          ', mean
  write ( *, '(a,g14.6)' ) '  PDF parameter VARIANCE =      ', variance

  do i = 1, sample_num
    call triangle_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test145

!*****************************************************************************80
!
!! TEST145 tests TRIANGULAR_CDF, TRIANGULAR_CDF_INV, TRIANGULAR_PDF;
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
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical triangular_check
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST145'
  write ( *, '(a)' ) '  For the Triangular PDF:'
  write ( *, '(a)' ) '  TRIANGULAR_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  TRIANGULAR_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  TRIANGULAR_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =      ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =      ', b

  if ( .not. triangular_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST145 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call triangular_sample ( a, b, seed, x )

    call triangular_pdf ( x, a, b, pdf )

    call triangular_cdf ( x, a, b, cdf )

    call triangular_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test146

!*****************************************************************************80
!
!! TEST146 tests TRIANGULAR_MEAN, TRIANGULAR_SAMPLE, TRIANGULAR_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  logical triangular_check
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST146'
  write ( *, '(a)' ) '  For the Triangular PDF:'
  write ( *, '(a)' ) '  TRIANGULAR_MEAN computes mean;'
  write ( *, '(a)' ) '  TRIANGULAR_SAMPLE samples;'
  write ( *, '(a)' ) '  TRIANGULAR_VARIANCE computes variance.'

  a = 1.0D+00
  b = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. triangular_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST146 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call triangular_mean ( a, b, mean )
  call triangular_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call triangular_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test147

!*****************************************************************************80
!
!! TEST147 tests UNIFORM_01_ORDER_SAMPLE;
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

  integer ( kind = 4 ), parameter :: n = 10

  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST147'
  write ( *, '(a)' ) '  For the Uniform 01 Order PDF:'
  write ( *, '(a)') '  UNIFORM_ORDER_SAMPLE samples.'
  write ( *, '(a)' ) ' '

  call uniform_01_order_sample ( n, seed, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Ordered sample:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,g14.6)' ) i, x(i)
  end do

  return
end
subroutine test148

!*****************************************************************************80
!
!! TEST148 tests UNIFORM_NSPHERE_SAMPLE;
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

  integer ( kind = 4 ), parameter :: n = 3

  integer ( kind = 4 ) i
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST148'
  write ( *, '(a)' ) '  For the Uniform PDF on the N-Sphere:'
  write ( *, '(a)' ) '  UNIFORM_NSPHERE_SAMPLE samples.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Dimension N of sphere =       ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Points on the sphere:'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    call uniform_nsphere_sample ( n, seed, x )
    write ( *, '(2x,i8,3g14.6)' ) i, x(1:n)
  end do

  return
end
subroutine test1485

!*****************************************************************************80
!
!! TEST1485 tests UNIFORM_01_CDF, UNIFORM_01_CDF_INV, UNIFORM_01_PDF;
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

  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) uniform_01_sample
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1485'
  write ( *, '(a)' ) '  For the Uniform 01 PDF:'
  write ( *, '(a)' ) '  UNIFORM_01_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  UNIFORM_01_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  UNIFORM_01_PDF evaluates the PDF;'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    x = uniform_01_sample ( seed )

    call uniform_01_pdf ( x, pdf )

    call uniform_01_cdf ( x, cdf )

    call uniform_01_cdf_inv ( cdf, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test1486

!*****************************************************************************80
!
!! TEST1486 tests UNIFORM_01_MEAN, UNIFORM_01_SAMPLE, UNIFORM_01_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) uniform_01_sample
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1486'
  write ( *, '(a)' ) '  For the Uniform 01 PDF:'
  write ( *, '(a)' ) '  UNIFORM_01_MEAN computes mean;'
  write ( *, '(a)' ) '  UNIFORM_01_SAMPLE samples;'
  write ( *, '(a)' ) '  UNIFORM_01_VARIANCE computes variance.'

  call uniform_01_mean ( mean )
  call uniform_01_variance ( variance )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF mean =            ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =        ', variance

  do i = 1, sample_num
    x(i) = uniform_01_sample ( seed )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test149

!*****************************************************************************80
!
!! TEST149 tests UNIFORM_CDF, UNIFORM_CDF_INV, UNIFORM_PDF;
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical uniform_check
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST149'
  write ( *, '(a)' ) '  For the Uniform PDF:'
  write ( *, '(a)' ) '  UNIFORM_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  UNIFORM_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  UNIFORM_PDF evaluates the PDF;'

  a = 1.0D+00
  b = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =      ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =      ', b

  if ( .not. uniform_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST149 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call uniform_sample ( a, b, seed, x )

    call uniform_pdf ( x, a, b, pdf )

    call uniform_cdf ( x, a, b, cdf )

    call uniform_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test150

!*****************************************************************************80
!
!! TEST150 tests UNIFORM_MEAN, UNIFORM_SAMPLE, UNIFORM_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  logical uniform_check
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST150'
  write ( *, '(a)' ) '  For the Uniform PDF:'
  write ( *, '(a)' ) '  UNIFORM_MEAN computes mean;'
  write ( *, '(a)' ) '  UNIFORM_SAMPLE samples;'
  write ( *, '(a)' ) '  UNIFORM_VARIANCE computes variance.'

  a = 1.0D+00
  b = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. uniform_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST150 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call uniform_mean ( a, b, mean )
  call uniform_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =        ', variance

  do i = 1, sample_num
    call uniform_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test151

!*****************************************************************************80
!
!! TEST151 tests UNIFORM_DISCRETE_CDF, UNIFORM_DISCRETE_CDF_INV, UNIFORM_DISCRETE_PDF;
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

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical uniform_discrete_check
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST151'
  write ( *, '(a)' ) '  For the Uniform Discrete PDF:'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_PDF evaluates the PDF;'

  a = 1
  b = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A =             ', a
  write ( *, '(a,i8)'    ) '  PDF parameter B =             ', b

  if ( .not. uniform_discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST151 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call uniform_discrete_sample ( a, b, seed, x )

    call uniform_discrete_pdf ( x, a, b, pdf )

    call uniform_discrete_cdf ( x, a, b, cdf )

    call uniform_discrete_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test152

!*****************************************************************************80
!
!! TEST152 tests UNIFORM_DISCRETE_MEAN, UNIFORM_DISCRETE_SAMPLE, UNIFORM_DISCRETE_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  logical uniform_discrete_check
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST152'
  write ( *, '(a)' ) '  For the Uniform discrete PDF:'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_MEAN computes the mean;'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_SAMPLE samples;'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_VARIANCE computes the variance.'

  a = 1
  b = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  PDF parameter A =             ', a
  write ( *, '(a,i8)'    ) '  PDF parameter B =             ', b

  if ( .not. uniform_discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST143 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call uniform_discrete_mean ( a, b, mean )
  call uniform_discrete_variance ( a, b, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call uniform_discrete_sample ( a, b, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test153

!*****************************************************************************80
!
!! TEST153 tests UNIFORM_DISCRETE_CDF, UNIFORM_DISCRETE_PDF.
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

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  logical uniform_discrete_check
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST153'
  write ( *, '(a)' ) '  For the Uniform discrete PDF.'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  UNIFORM_DISCRETE_CDF evaluates the CDF.'

  a = 1
  b = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  PDF parameter A =             ', a
  write ( *, '(a,i8)' ) '  PDF parameter B =             ', b

  if ( .not. uniform_discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST153 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      PDF(X)      CDF(X)'
  write ( *, '(a)' ) ' '

  do x = 0, 6
    call uniform_discrete_pdf ( x, a, b, pdf )
    call uniform_discrete_cdf ( x, a, b, cdf )
    write ( *, '(2x,i8,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test154

!*****************************************************************************80
!
!! TEST154 tests VON_MISES_CDF, VON_MISES_CDF_INV, VON_MISES_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical von_mises_check
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST154'
  write ( *, '(a)' ) '  For the Von Mises PDF:'
  write ( *, '(a)' ) '  VON_MISES_CDF evaluates the CDF.'
  write ( *, '(a)' ) '  VON_MISES_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  VON_MISES_PDF evaluates the PDF.'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =      ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =      ', b

  if ( .not. von_mises_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST154 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call von_mises_sample ( a, b, seed, x )

    call von_mises_pdf ( x, a, b, pdf )

    call von_mises_cdf ( x, a, b, cdf )

    call von_mises_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test155

!*****************************************************************************80
!
!! TEST155 tests VON_MISES_MEAN, VON_MISES_SAMPLE, VON_MISES_CIRCULAR_VARIANCE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) circular_variance
  logical von_mises_check
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST155'
  write ( *, '(a)' ) '  For the Von Mises PDF:'
  write ( *, '(a)' ) '  VON_MISES_MEAN computes the mean;'
  write ( *, '(a)' ) '  VON_MISES_SAMPLE samples.'
  write ( *, '(a)' ) &
    '  VON_MISES_CIRCULAR_VARIANCE computes the circular_variance.'

  a = 1.0D+00
  b = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. von_mises_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST155 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call von_mises_mean ( a, b, mean )
  call von_mises_circular_variance ( a, b, circular_variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF circular variance =       ', circular_variance

  do i = 1, sample_num
    call von_mises_sample ( a, b, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_circular_variance ( sample_num, x, circular_variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =              ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =              ', mean
  write ( *, '(a,g14.6)' ) '  Sample circular variance = ', circular_variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =           ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =           ', xmin

  return
end
subroutine test1555

!*****************************************************************************80
!
!! TEST1555 tests VON_MISES_CDF, VON_MISES_CDF_VALUES.
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
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST1555:'
  write ( *, '(a)' ) '  VON_MISES_CDF evaluates the cumulative distribution'
  write ( *, '(a)' ) '    function for the von Mises PDF.'
  write ( *, '(a)' ) '  VON_MISES_CDF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A is the dominant angle;'
  write ( *, '(a)' ) '  B is a measure of spread;'
  write ( *, '(a)' ) '  X is the angle;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A         B         X          Exact F                Computed F'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call von_mises_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call von_mises_cdf ( x, a, b, fx2 )

    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4,g24.16,g24.16)' ) a, b, x, fx, fx2

  end do

  return
end
subroutine test156

!*****************************************************************************80
!
!! TEST156 tests WEIBULL_CDF, WEIBULL_CDF_INV, WEIBULL_PDF.
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
  real ( kind = 8 ) c
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical weibull_check
  real ( kind = 8 ) x
  real ( kind = 8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST156'
  write ( *, '(a)' ) '  For the Weibull PDF:'
  write ( *, '(a)' ) '  WEIBULL_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  WEIBULL_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  WEIBULL_PDF evaluates the PDF;'

  x = 3.0D+00

  a = 2.0D+00
  b = 3.0D+00
  c = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. weibull_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST156 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call weibull_sample ( a, b, c, seed, x )

    call weibull_pdf ( x, a, b, c, pdf )

    call weibull_cdf ( x, a, b, c, cdf )

    call weibull_cdf_inv ( cdf, a, b, c, x2 )

    write ( *, '(2x,4g14.6)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test157

!*****************************************************************************80
!
!! TEST157 tests WEIBULL_MEAN, WEIBULL_SAMPLE, WEIBULL_VARIANCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  logical weibull_check
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST157'
  write ( *, '(a)' ) '  For the Weibull PDF:'
  write ( *, '(a)' ) '  WEIBULL_MEAN computes the mean;'
  write ( *, '(a)' ) '  WEIBULL_SAMPLE samples;'
  write ( *, '(a)' ) '  WEIBULL_VARIANCE computes the variance.'

  a = 2.0D+00
  b = 3.0D+00
  c = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =       ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b
  write ( *, '(a,g14.6)' ) '  PDF parameter C =             ', c

  if ( .not. weibull_check ( a, b, c ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST157 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call weibull_mean ( a, b, c, mean )
  call weibull_variance ( a, b, c, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call weibull_sample ( a, b, c, seed, x(i) )
  end do

  call r8vec_mean ( sample_num, x, mean )
  call r8vec_variance ( sample_num, x, variance )
  call r8vec_max ( sample_num, x, imax, xmax )
  call r8vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine test158

!*****************************************************************************80
!
!! TEST158 tests WEIBULL_DISCRETE_CDF, WEIBULL_DISCRETE_CDF_INV, WEIBULL_DISCRETE_PDF.
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
  real ( kind = 8 ) cdf
  integer ( kind = 4 ) i
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) :: seed = 123456789
  logical weibull_discrete_check
  integer ( kind = 4 ) x
  integer ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST158'
  write ( *, '(a)' ) '  For the Weibull Discrete PDF,'
  write ( *, '(a)' ) '  WEIBULL_DISCRETE_CDF evaluates the CDF;'
  write ( *, '(a)' ) '  WEIBULL_DISCRETE_CDF_INV inverts the CDF.'
  write ( *, '(a)' ) '  WEIBULL_DISCRETE_PDF evaluates the PDF;'

  a = 0.50D+00
  b = 1.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =             ', b

  if ( .not. weibull_discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST158 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X            PDF           CDF            CDF_INV'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    call weibull_discrete_sample ( a, b, seed, x )

    call weibull_discrete_pdf ( x, a, b, pdf )

    call weibull_discrete_cdf ( x, a, b, cdf )

    call weibull_discrete_cdf_inv ( cdf, a, b, x2 )

    write ( *, '(2x,i14,2g14.6,i14)' ) x, pdf, cdf, x2

  end do

  return
end
subroutine test159

!*****************************************************************************80
!
!! TEST159 tests WEIBULL_DISCRETE_CDF, WEIBULL_DISCRETE_PDF.
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
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  logical weibull_discrete_check
  integer ( kind = 4 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST159'
  write ( *, '(a)' ) '  For the Weibull Discrete PDF:'
  write ( *, '(a)' ) '  WEIBULL_DISCRETE_PDF evaluates the PDF;'
  write ( *, '(a)' ) '  WEIBULL_DISCRETE_CDF evaluates the CDF.'

  a = 0.50D+00
  b = 1.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =     ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =     ', b

  if ( .not. weibull_discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST159 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      PDF(X)      CDF(X)'
  write ( *, '(a)' ) ' '

  do x = 0, 10
    call weibull_discrete_pdf ( x, a, b, pdf )
    call weibull_discrete_cdf ( x, a, b, cdf )
    write ( *, '(2x,i8,2g14.6)' ) x, pdf, cdf
  end do

  return
end
subroutine test160

!*****************************************************************************80
!
!! TEST160 tests WEIBULL_DISCRETE_SAMPLE.
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

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  logical weibull_discrete_check
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST160'
  write ( *, '(a)' ) '  For the discrete Weibull PDF:'
  write ( *, '(a)' ) '  WEIBULL_DISCRETE_SAMPLE samples.'

  a = 0.5D+00
  b = 1.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =     ', a
  write ( *, '(a,g14.6)' ) '  PDF parameter B =     ', b

  if ( .not. weibull_discrete_check ( a, b ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST160 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  do i = 1, sample_num
    call weibull_discrete_sample ( a, b, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
subroutine test161

!*****************************************************************************80
!
!! TEST161 tests ZIPF_CDF and ZIPF_PDF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) cdf
  real ( kind = 8 ) pdf
  integer ( kind = 4 ) x
  logical zipf_check

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST161'
  write ( *, '(a)' ) '  For the Zipf PDF:'
  write ( *, '(a)' ) '  ZIPF_PDF evaluates the PDF.'
  write ( *, '(a)' ) '  ZIPF_CDF evaluates the CDF.'

  a = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A = ', a

  if ( .not. zipf_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST161 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  X    PDF(X)       CDF(X)'
  write ( *, '(a)' ) ' '

  do x = 1, 20

    call zipf_pdf ( x, a, pdf )
    call zipf_cdf ( x, a, cdf )
    write ( *, '(2x,i8,2x,2g14.6)' ) x, pdf, cdf

  end do

  return
end
subroutine test162

!*****************************************************************************80
!
!! TEST162 tests ZIPF_SAMPLE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  real ( kind = 8 ) mean
  integer ( kind = 4 ) :: seed = 123456789
  real ( kind = 8 ) variance
  integer ( kind = 4 ) x(sample_num)
  integer ( kind = 4 ) xmax
  integer ( kind = 4 ) xmin
  logical zipf_check

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST162'
  write ( *, '(a)' ) '  For the Zipf PDF:'
  write ( *, '(a)' ) '  ZIPF_SAMPLE samples.'

  a = 4.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  PDF parameter A =             ', a

  if ( .not. zipf_check ( a ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TEST162 - Fatal error!'
    write ( *, '(a)' ) '  The parameters are not legal.'
    return
  end if

  call zipf_mean ( a, mean )
  call zipf_variance ( a, variance )

  write ( *, '(a,g14.6)' ) '  PDF mean =                    ', mean
  write ( *, '(a,g14.6)' ) '  PDF variance =                ', variance

  do i = 1, sample_num
    call zipf_sample ( a, seed, x(i) )
  end do

  call i4vec_mean ( sample_num, x, mean )
  call i4vec_variance ( sample_num, x, variance )
  call i4vec_max ( sample_num, x, imax, xmax )
  call i4vec_min ( sample_num, x, imin, xmin )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,i8)'    ) '  Sample maximum =  ', xmax
  write ( *, '(a,i8)'    ) '  Sample minimum =  ', xmin

  return
end
