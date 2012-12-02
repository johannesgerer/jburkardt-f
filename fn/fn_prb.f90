program main

!*****************************************************************************80
!
!! MAIN is the main program for FN_PRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FN_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( * ,'(a)' ) '  Test the FN library.'

  call acos_test ( )
  call acosh_test ( )
  call ai_test ( )
  call aid_test ( )
  call asin_test ( )
  call asinh_test ( )
  call atan_test ( )
  call atan2_test ( )
  call atanh_test ( )
  call besi0_test ( )
  call besi1_test ( )
  call besj0_test ( )
  call besj1_test ( )
  call besk_test ( )
  call besk0_test ( )
  call besk1_test ( )
  call besy0_test ( )
  call besy1_test ( )
  call beta_test ( )
  call betai_test ( )
  call bi_test ( )
  call bid_test ( )
  call binom_test ( )
  call cbrt_test ( )
  call chi_test ( )
  call chu_test ( )
  call ci_test ( )
  call cin_test ( )
  call cinh_test ( )
  call cos_test ( )
  call cos_deg_test ( )
  call cosh_test ( )
  call cot_test ( )
  call dawson_test ( )
  call e1_test ( )
  call ei_test ( )
  call erf_test ( )
  call erfc_test ( )
  call exp_test ( )
  call fac_test ( )
  call gamma_test ( )
  call gamma_inc_test ( )
  call gamma_inc_tricomi_test ( )
  call lbeta_test ( )
  call li_test ( )
  call lngam_test ( )
  call log_test ( )
  call log10_test ( )
  call poch_test ( )
  call psi_test ( )
  call rand_test ( )
  call shi_test ( )
  call si_test ( )
  call sin_test ( )
  call sin_deg_test ( )
  call sinh_test ( )
  call spence_test ( )
  call sqrt_test ( )
  call tan_test ( )
  call tanh_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FN_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine acos_test ( )

!*****************************************************************************80
!
!! ACOS_TEST tests R4_ACOS and R8_ACOS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_acos
  real ( kind = 8 ) r8_acos
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ACOS_TEST:'
  write ( *, '(a)' ) '  Test ARCCOS_VALUES, R4_ACOS, R8_ACOS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X      ARCCOS(X)'
  write ( *, '(a)' ) '                   R4_ACOS(X)         Diff'
  write ( *, '(a)' ) '                   R8_ACOS(X)         Diff'

  n_data = 0

  do

    call arccos_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_acos ( real ( x, kind = 4 ) )
    fx3 = r8_acos ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine acosh_test ( )

!*****************************************************************************80
!
!! ACOSH_TEST tests R4_ACOSH and R8_ACOSH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_acosh
  real ( kind = 8 ) r8_acosh
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ACOSH_TEST:'
  write ( *, '(a)' ) '  Test ARCCOSH_VALUES, R4_ACOSH, R8_ACOSH'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X      ARCCOSH(X)'
  write ( *, '(a)' ) '                   R4_ACOSH(X)        Diff'
  write ( *, '(a)' ) '                   R8_ACOSH(X)        Diff'

  n_data = 0

  do

    call arccosh_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_acosh ( real ( x, kind = 4 ) )
    fx3 = r8_acosh ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine ai_test ( )

!*****************************************************************************80
!
!! AI_TEST tests R4_AI and R8_AI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_ai
  real ( kind = 8 ) r8_ai
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AI_TEST:'
  write ( *, '(a)' ) '  Test AIRY_AI_VALUES, R4_AI, R8_AI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X      AIRY_AI(X)'
  write ( *, '(a)' ) '                      R4_AI(X)        Diff'
  write ( *, '(a)' ) '                      R8_AI(X)        Diff'

  n_data = 0

  do

    call airy_ai_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_ai ( real ( x, kind = 4 ) )
    fx3 = r8_ai ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine aid_test ( )

!*****************************************************************************80
!
!! AID_TEST tests R4_AID and R8_AID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_aid
  real ( kind = 8 ) r8_aid
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AID_TEST:'
  write ( *, '(a)' ) '  Test AIRY_AI_PRIME_VALUES, R4_AID, R8_AID'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X     AIRY_AID(X)'
  write ( *, '(a)' ) '                     R4_AID(X)        Diff'
  write ( *, '(a)' ) '                     R8_AID(X)        Diff'

  n_data = 0

  do

    call airy_ai_prime_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_aid ( real ( x, kind = 4 ) )
    fx3 = r8_aid ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine asin_test ( )

!*****************************************************************************80
!
!! ASIN_TEST tests R4_ASIN and R8_ASIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_asin
  real ( kind = 8 ) r8_asin
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASIN_TEST:'
  write ( *, '(a)' ) '  Test ARCSIN_VALUES, R4_ASIN, R8_ASIN'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X      ARCSIN(X)'
  write ( *, '(a)' ) '                   R4_ASIN(X)         Diff'
  write ( *, '(a)' ) '                   R8_ASIN(X)         Diff'

  n_data = 0

  do

    call arcsin_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_asin ( real ( x, kind = 4 ) )
    fx3 = r8_asin ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine asinh_test ( )

!*****************************************************************************80
!
!! ASINH_TEST tests R4_ASINH and R8_ASINH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_asinh
  real ( kind = 8 ) r8_asinh
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASINH_TEST:'
  write ( *, '(a)' ) '  Test ARCSINH_VALUES, R4_ASINH, R8_ASINH'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X     ARCSINH(X)'
  write ( *, '(a)' ) '                  R4_ASINH(X)         Diff'
  write ( *, '(a)' ) '                  R8_ASINH(X)         Diff'

  n_data = 0

  do

    call arcsinh_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_asinh ( real ( x, kind = 4 ) )
    fx3 = r8_asinh ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine atan_test ( )

!*****************************************************************************80
!
!! ATAN_TEST tests R4_ATAN and R8_ATAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_atan
  real ( kind = 8 ) r8_atan
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ATAN_TEST:'
  write ( *, '(a)' ) '  Test ARCTAN_VALUES, R4_ATAN, R8_ATAN'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X      ARCTAN(X)'
  write ( *, '(a)' ) '                   R4_ATAN(X)         Diff'
  write ( *, '(a)' ) '                   R8_ATAN(X)         Diff'

  n_data = 0

  do

    call arctan_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_atan ( real ( x, kind = 4 ) )
    fx3 = r8_atan ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine atan2_test ( )

!*****************************************************************************80
!
!! ATAN2_TEST tests R4_ATAN2 and R8_ATAN2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_atan2
  real ( kind = 8 ) r8_atan2
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ATAN2_TEST:'
  write ( *, '(a)' ) '  Test ARCTAN2_VALUES, R4_ATAN2, R8_ATAN2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X               Y   ARCTAN2(Y,X)'
  write ( *, '(a)' ) '                                R4_ATAN2(Y,X)         Diff'
  write ( *, '(a)' ) '                                R8_ATAN2(Y,X)         Diff'

  n_data = 0

  do

    call arctan2_values ( n_data, x, y, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_atan2 ( real ( y, kind = 4 ), real ( x, kind = 4 ) )
    fx3 = r8_atan2 ( y, x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,f14.4,2x,g14.6)' )     x, y, fx1
    write ( *, '(32x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(32x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine atanh_test ( )

!*****************************************************************************80
!
!! ATANH_TEST tests R4_ATANH and R8_ATANH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_atanh
  real ( kind = 8 ) r8_atanh
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ATANH_TEST:'
  write ( *, '(a)' ) '  Test ARCTANH_VALUES, R4_ATANH, R8_ATANH'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X     ARCTANH(X)'
  write ( *, '(a)' ) '                  R4_ATANH(X)         Diff'
  write ( *, '(a)' ) '                  R8_ATANH(X)         Diff'

  n_data = 0

  do

    call arctanh_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_atanh ( real ( x, kind = 4 ) )
    fx3 = r8_atanh ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besi0_test ( )

!*****************************************************************************80
!
!! BESI0_TEST tests R4_BESI0 and R8_BESI0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besi0
  real ( kind = 8 ) r8_besi0
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESI0_TEST:'
  write ( *, '(a)' ) '  Test BESI0_VALUES, R4_BESI0, R8_BESI0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESI0(X)'
  write ( *, '(a)' ) '                  R4_BESI0(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESI0(X)         Diff'

  n_data = 0

  do

    call bessel_i0_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besi0 ( real ( x, kind = 4 ) )
    fx3 = r8_besi0 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besi1_test ( )

!*****************************************************************************80
!
!! BESI1_TEST tests R4_BESI1 and R8_BESI1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besi1
  real ( kind = 8 ) r8_besi1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESI1_TEST:'
  write ( *, '(a)' ) '  Test BESI1_VALUES, R4_BESI1, R8_BESI1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESI1(X)'
  write ( *, '(a)' ) '                  R4_BESI1(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESI1(X)         Diff'

  n_data = 0

  do

    call bessel_i1_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besi1 ( real ( x, kind = 4 ) )
    fx3 = r8_besi1 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besj0_test ( )

!*****************************************************************************80
!
!! BESJ0_TEST tests R4_BESJ0 and R8_BESJ0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besj0
  real ( kind = 8 ) r8_besj0
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESJ0_TEST:'
  write ( *, '(a)' ) '  Test BESJ0_VALUES, R4_BESJ0, R8_BESJ0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESJ0(X)'
  write ( *, '(a)' ) '                  R4_BESJ0(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESJ0(X)         Diff'

  n_data = 0

  do

    call bessel_j0_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besj0 ( real ( x, kind = 4 ) )
    fx3 = r8_besj0 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besj1_test ( )

!*****************************************************************************80
!
!! BESJ1_TEST tests R4_BESJ1 and R8_BESJ1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besj1
  real ( kind = 8 ) r8_besj1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESJ1_TEST:'
  write ( *, '(a)' ) '  Test BESJ1_VALUES, R4_BESJ1, R8_BESJ1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESJ1(X)'
  write ( *, '(a)' ) '                  R4_BESJ1(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESJ1(X)         Diff'

  n_data = 0

  do

    call bessel_j1_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besj1 ( real ( x, kind = 4 ) )
    fx3 = r8_besj1 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besk_test ( )

!*****************************************************************************80
!
!! BESK_TEST tests R4_BESK and R8_BESK.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 November 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 4 ) r4_besk
  real ( kind = 8 ) r8_besk
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESK_TEST:'
  write ( *, '(a)' ) '  Test BESK_KX_VALUES, R4_BESK, R8_BESK'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '              NU             X       BESK(X)'
  write ( *, '(a)' ) '                                  R4_BESK(X)         Diff'
  write ( *, '(a)' ) '                                  R8_BESK(X)         Diff'

  n_data = 0

  do

    call bessel_kx_values ( n_data, nu, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besk ( real ( nu, kind = 4 ), real ( x, kind = 4 ) )
    fx3 = r8_besk ( nu, x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,f14.4,2x,g14.6)' )     nu, x, fx1
    write ( *, '(16x,16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besk0_test ( )

!*****************************************************************************80
!
!! BESK0_TEST tests R4_BESK0 and R8_BESK0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besk0
  real ( kind = 8 ) r8_besk0
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESK0_TEST:'
  write ( *, '(a)' ) '  Test BESK0_VALUES, R4_BESK0, R8_BESK0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESK0(X)'
  write ( *, '(a)' ) '                  R4_BESK0(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESK0(X)         Diff'

  n_data = 0

  do

    call bessel_k0_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besk0 ( real ( x, kind = 4 ) )
    fx3 = r8_besk0 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besk1_test ( )

!*****************************************************************************80
!
!! BESK1_TEST tests R4_BESK1 and R8_BESK1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besk1
  real ( kind = 8 ) r8_besk1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESK1_TEST:'
  write ( *, '(a)' ) '  Test BESK1_VALUES, R4_BESK1, R8_BESK1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESK1(X)'
  write ( *, '(a)' ) '                  R4_BESK1(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESK1(X)         Diff'

  n_data = 0

  do

    call bessel_k1_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besk1 ( real ( x, kind = 4 ) )
    fx3 = r8_besk1 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besy0_test ( )

!*****************************************************************************80
!
!! BESY0_TEST tests R4_BESY0 and R8_BESY0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besy0
  real ( kind = 8 ) r8_besy0
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESY0_TEST:'
  write ( *, '(a)' ) '  Test BESY0_VALUES, R4_BESY0, R8_BESY0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESY0(X)'
  write ( *, '(a)' ) '                  R4_BESY0(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESY0(X)         Diff'

  n_data = 0

  do

    call bessel_y0_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besy0 ( real ( x, kind = 4 ) )
    fx3 = r8_besy0 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine besy1_test ( )

!*****************************************************************************80
!
!! BESY1_TEST tests R4_BESY1 and R8_BESY1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_besy1
  real ( kind = 8 ) r8_besy1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESY1_TEST:'
  write ( *, '(a)' ) '  Test BESY1_VALUES, R4_BESY1, R8_BESY1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       BESY1(X)'
  write ( *, '(a)' ) '                  R4_BESY1(X)         Diff'
  write ( *, '(a)' ) '                  R8_BESY1(X)         Diff'

  n_data = 0

  do

    call bessel_y1_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_besy1 ( real ( x, kind = 4 ) )
    fx3 = r8_besy1 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine beta_test ( )

!*****************************************************************************80
!
!! BETA_TEST tests R4_BETA and R8_BETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_beta
  real ( kind = 8 ) r8_beta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_TEST:'
  write ( *, '(a)' ) '  Test BETA_VALUES, R4_BETA, R8_BETA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             A               B      BETA(A,B)'
  write ( *, '(a)' ) '                                 R4_BETA(A,B)         Diff'
  write ( *, '(a)' ) '                                 R8_BETA(A,B)         Diff'

  n_data = 0

  do

    call beta_values ( n_data, a, b, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, b

    fx2 = r4_beta ( real ( a, kind = 4 ), real ( b, kind = 4 ) )
    fx3 = r8_beta ( a, b )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, b, fx1
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( fx1 - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine betai_test ( )

!*****************************************************************************80
!
!! BETAI_TEST tests R4_BETAI and R8_BETAI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_betai
  real ( kind = 8 ) r8_betai
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETAI_TEST:'
  write ( *, '(a)' ) '  Test BETA_INC_VALUES, R4_BETAI, R8_BETAI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X               A               B     BETAI(A,B)'
  write ( *, '(a)' ) &
    '                                                R4_BETAI(A,B)' &
    // '         Diff'
  write ( *, '(a)' ) &
    '                                                R8_BETAI(A,B)' &
    // '         Diff'

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_betai ( real ( x, kind = 4 ), real ( a ), real ( b ) )
    fx3 = r8_betai ( x, a, b )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,f14.4,2x,g14.4,2x,g14.6)' ) x, a, b, fx1
    write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( fx1 - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine bi_test ( )

!*****************************************************************************80
!
!! BI_TEST tests R4_BI and R8_BI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_bi
  real ( kind = 8 ) r8_bi
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BI_TEST:'
  write ( *, '(a)' ) '  Test AIRY_BI_VALUES, R4_BI, R8_BI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X      AIRY_BI(X)'
  write ( *, '(a)' ) '                      R4_BI(X)        Diff'
  write ( *, '(a)' ) '                      R8_BI(X)        Diff'

  n_data = 0

  do

    call airy_bi_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_bi ( real ( x, kind = 4 ) )
    fx3 = r8_bi ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine bid_test ( )

!*****************************************************************************80
!
!! BID_TEST tests R4_BID and R8_BID.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 April 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_bid
  real ( kind = 8 ) r8_bid
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BID_TEST:'
  write ( *, '(a)' ) '  Test AIRY_BI_PRIME_VALUES, R4_BID, R8_BID'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X     AIRY_BID(X)'
  write ( *, '(a)' ) '                     R4_BID(X)        Diff'
  write ( *, '(a)' ) '                     R8_BID(X)        Diff'

  n_data = 0

  do

    call airy_bi_prime_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_bid ( real ( x, kind = 4 ) )
    fx3 = r8_bid ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine binom_test ( )

!*****************************************************************************80
!
!! BINOM_TEST tests R4_BINOM and R8_BINOM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real ( kind = 8 ) diff
  integer ( kind = 4 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_binom
  real ( kind = 8 ) r8_binom

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BINOM_TEST:'
  write ( *, '(a)' ) '  Test BINOM_VALUES, R4_BINOM, R8_BINOM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             A               B     BINOM(A,B)'
  write ( *, '(a)' ) '                                R4_BINOM(A,B)         Diff'
  write ( *, '(a)' ) '                                R8_BINOM(A,B)         Diff'

  n_data = 0

  do

    call binomial_values ( n_data, a, b, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_binom ( a, b )
    fx3 = r8_binom ( a, b )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i14,2x,i14,2x,i14)' )  a, b, fx1
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( real ( fx1 ) - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( dble ( fx1 ) - fx3 )

  end do

  return
end
subroutine cbrt_test ( )

!*****************************************************************************80
!
!! CBRT_TEST tests R4_CBRT and R8_CBRT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_cbrt
  real ( kind = 8 ) r8_cbrt
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CBRT_TEST:'
  write ( *, '(a)' ) '  Test CBRT_VALUES, R4_CBRT, R8_CBRT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X        CBRT(X)'
  write ( *, '(a)' ) '                   R4_CBRT(X)         Diff'
  write ( *, '(a)' ) '                   R8_CBRT(X)         Diff'

  n_data = 0

  do

    call cbrt_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_cbrt ( real ( x, kind = 4 ) )
    fx3 = r8_cbrt ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,g14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine chu_test ( )

!*****************************************************************************80
!
!! CHU_TEST tests R4_CHU and R8_CHU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_chu
  real ( kind = 8 ) r8_chu
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHU_TEST:'
  write ( *, '(a)' ) '  Test HYPERGEOMETRIC_U_VALUES, R4_CHU, R8_CHU.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             A               B               X     CHU(A,B,X)'
  write ( *, '(a)' ) &
   '                                                R4_CHU(A,B,X)' &
    // '         Diff'
  write ( *, '(a)' ) &
    '                                                R8_CHU(A,B,X)' &
    // '         Diff'

  n_data = 0

  do

    call hypergeometric_u_values ( n_data, a, b, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_chu ( real ( a ), real ( b ), real ( x, kind = 4 ) )
    fx3 = r8_chu ( a, b, x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,f14.4,2x,g14.4,2x,g14.6)' ) a, b, x, fx1
    write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( fx1 - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine chi_test ( )

!*****************************************************************************80
!
!! CHI_TEST tests R4_CHI and R8_CHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_chi
  real ( kind = 8 ) r8_chi
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHI_TEST:'
  write ( *, '(a)' ) '  Test CHI_VALUES, R4_CHI, R8_CHI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X          CHI(X)'
  write ( *, '(a)' ) '                     R4_CHI(X)        Diff'
  write ( *, '(a)' ) '                     R8_CHI(X)        Diff'

  n_data = 0

  do

    call chi_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_chi ( real ( x, kind = 4 ) )
    fx3 = r8_chi ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine ci_test ( )

!*****************************************************************************80
!
!! CI_TEST tests R4_CI and R8_CI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_ci
  real ( kind = 8 ) r8_ci
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CI_TEST:'
  write ( *, '(a)' ) '  Test CI_VALUES, R4_CI, R8_CI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X           CI(X)'
  write ( *, '(a)' ) '                      R4_CI(X)        Diff'
  write ( *, '(a)' ) '                      R8_CI(X)        Diff'

  n_data = 0

  do

    call ci_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_ci ( real ( x, kind = 4 ) )
    fx3 = r8_ci ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine cin_test ( )

!*****************************************************************************80
!
!! CIN_TEST tests R4_CIN and R8_CIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_cin
  real ( kind = 8 ) r8_cin
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIN_TEST:'
  write ( *, '(a)' ) '  Test CIN_VALUES, R4_CIN, R8_CIN.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X          CIN(X)'
  write ( *, '(a)' ) '                     R4_CIN(X)        Diff'
  write ( *, '(a)' ) '                     R8_CIN(X)        Diff'

  n_data = 0

  do

    call cin_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_cin ( real ( x, kind = 4 ) )
    fx3 = r8_cin ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine cinh_test ( )

!*****************************************************************************80
!
!! CINH_TEST tests R4_CINH and R8_CINH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_cinh
  real ( kind = 8 ) r8_cinh
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CINH_TEST:'
  write ( *, '(a)' ) '  Test CINH_VALUES, R4_CINH, R8_CINH.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         CINH(X)'
  write ( *, '(a)' ) '                    R4_CINH(X)        Diff'
  write ( *, '(a)' ) '                    R8_CINH(X)        Diff'

  n_data = 0

  do

    call cinh_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_cinh ( real ( x, kind = 4 ) )
    fx3 = r8_cinh ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine cos_test ( )

!*****************************************************************************80
!
!! COS_TEST tests R4_COS and R8_COS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_cos
  real ( kind = 8 ) r8_cos
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COS_TEST:'
  write ( *, '(a)' ) '  Test COS_VALUES, R4_COS, R8_COS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         COS(X)'
  write ( *, '(a)' ) '                    R4_COS(X)         Diff'
  write ( *, '(a)' ) '                    R8_COS(X)         Diff'

  n_data = 0

  do

    call cos_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_cos ( real ( x, kind = 4 ) )
    fx3 = r8_cos ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine cos_deg_test ( )

!*****************************************************************************80
!
!! COS_DEG_TEST tests R4_COS_DEG and R8_COS_DEG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_cos_deg
  real ( kind = 8 ) r8_cos_deg
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COS_DEG_TEST:'
  write ( *, '(a)' ) '  Test COS_DEGREE_VALUES, R4_COS_DEG, R8_COS_DEG.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X     COS_DEG(X)'
  write ( *, '(a)' ) '                R4_COS_DEG(X)         Diff'
  write ( *, '(a)' ) '                R8_COS_DEG(X)         Diff'

  n_data = 0

  do

    call cos_degree_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_cos_deg ( real ( x, kind = 4 ) )
    fx3 = r8_cos_deg ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine cosh_test ( )

!*****************************************************************************80
!
!! COSH_TEST tests R4_COSH and R8_COSH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_cosh
  real ( kind = 8 ) r8_cosh
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COSH_TEST:'
  write ( *, '(a)' ) '  Test COSH_VALUES, R4_COSH, R8_COSH.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         COSH(X)'
  write ( *, '(a)' ) '                    R4_COSH(X)        Diff'
  write ( *, '(a)' ) '                    R8_COSH(X)        Diff'

  n_data = 0

  do

    call cosh_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_cosh ( real ( x, kind = 4 ) )
    fx3 = r8_cosh ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine cot_test ( )

!*****************************************************************************80
!
!! COT_TEST tests R4_COT and R8_COT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_cot
  real ( kind = 8 ) r8_cot
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COT_TEST:'
  write ( *, '(a)' ) '  Test COT_VALUES, R4_COT, R8_COT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         COT(X)'
  write ( *, '(a)' ) '                    R4_COT(X)         Diff'
  write ( *, '(a)' ) '                    R8_COT(X)         Diff'

  n_data = 0

  do

    call cot_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_cot ( real ( x, kind = 4 ) )
    fx3 = r8_cot ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine dawson_test ( )

!*****************************************************************************80
!
!! DAWSON_TEST tests R4_DAWSON and R8_DAWSON.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_dawson
  real ( kind = 8 ) r8_dawson
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DAWSON_TEST:'
  write ( *, '(a)' ) '  Test DAWSON_VALUES, R4_DAWSON, R8_DAWSON.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X      DAWSON(X)'
  write ( *, '(a)' ) '                 R4_DAWSON(X)         Diff'
  write ( *, '(a)' ) '                 R8_DAWSON(X)         Diff'

  n_data = 0

  do

    call dawson_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_dawson ( real ( x, kind = 4 ) )
    fx3 = r8_dawson ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine e1_test ( )

!*****************************************************************************80
!
!! E1_TEST tests R4_E1 and R8_E1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_e1
  real ( kind = 8 ) r8_e1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'E1_TEST:'
  write ( *, '(a)' ) '  Test E1_VALUES, R4_E1, R8_E1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X           E1(X)'
  write ( *, '(a)' ) '                      R4_E1(X)        Diff'
  write ( *, '(a)' ) '                      R8_E1(X)        Diff'

  n_data = 0

  do

    call e1_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_e1 ( real ( x, kind = 4 ) )
    fx3 = r8_e1 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine ei_test ( )

!*****************************************************************************80
!
!! EI_TEST tests R4_EI and R8_EI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_ei
  real ( kind = 8 ) r8_ei
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EI_TEST:'
  write ( *, '(a)' ) '  Test EI_VALUES, R4_EI, R8_EI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X           EI(X)'
  write ( *, '(a)' ) '                      R4_EI(X)        Diff'
  write ( *, '(a)' ) '                      R8_EI(X)        Diff'

  n_data = 0

  do

    call ei_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_ei ( real ( x, kind = 4 ) )
    fx3 = r8_ei ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine erf_test ( )

!*****************************************************************************80
!
!! ERF_TEST tests R4_ERF and R8_ERF.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_erf
  real ( kind = 8 ) r8_erf
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ERF_TEST:'
  write ( *, '(a)' ) '  Test ERF_VALUES, R4_ERF R8_ERF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X          ERF(X)'
  write ( *, '(a)' ) '                     R4_ERF(X)        Diff'
  write ( *, '(a)' ) '                     R8_ERF(X)        Diff'

  n_data = 0

  do

    call erf_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_erf ( real ( x, kind = 4 ) )
    fx3 = r8_erf ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine erfc_test ( )

!*****************************************************************************80
!
!! ERFC_TEST tests R4_ERFC and R8_ERFC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_erfc
  real ( kind = 8 ) r8_erfc
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ERFC_TEST:'
  write ( *, '(a)' ) '  Test ERFC_VALUES, R4_ERFC R8_ERFC'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         ERFC(X)'
  write ( *, '(a)' ) '                    R4_ERFC(X)        Diff'
  write ( *, '(a)' ) '                    R8_ERFC(X)        Diff'

  n_data = 0

  do

    call erfc_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_erfc ( real ( x, kind = 4 ) )
    fx3 = r8_erfc ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine exp_test ( )

!*****************************************************************************80
!
!! EXP_TEST tests R4_EXP and R8_EXP.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_exp
  real ( kind = 8 ) r8_exp
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXP_TEST:'
  write ( *, '(a)' ) '  Test EXP_VALUES, R4_EXP, R8_EXP'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X          EXP(X)'
  write ( *, '(a)' ) '                     R4_EXP(X)        Diff'
  write ( *, '(a)' ) '                     R8_EXP(X)        Diff'

  n_data = 0

  do

    call exp_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_exp ( real ( x, kind = 4 ) )
    fx3 = r8_exp ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine fac_test ( )

!*****************************************************************************80
!
!! FAC_TEST tests R4_FAC and R8_FAC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  integer ( kind = 4 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_fac
  real ( kind = 8 ) r8_fac

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FAC_TEST:'
  write ( *, '(a)' ) '  Test FACTORIAL_VALUES, R4_FAC, R8_FAC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             N          FAC(N)'
  write ( *, '(a)' ) '                     R4_FAC(N)        Diff'
  write ( *, '(a)' ) '                     R8_FAC(N)        Diff'

  n_data = 0

  do

    call factorial_values ( n_data, n, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_fac ( n )
    fx3 = r8_fac ( n )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i14,2x,i14)' )       n, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine gamma_test ( )

!*****************************************************************************80
!
!! GAMMA_TEST tests R4_GAMMA and R8_GAMMA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_gamma
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_TEST:'
  write ( *, '(a)' ) '  Test GAMMA_VALUES, R4_GAMMA, R8_GAMMA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X        GAMMA(X)'
  write ( *, '(a)' ) '                   R4_GAMMA(X)        Diff'
  write ( *, '(a)' ) '                   R8_GAMMA(X)        Diff'

  n_data = 0

  do

    call gamma_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_gamma ( real ( x, kind = 4 ) )
    fx3 = r8_gamma ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine gamma_inc_test ( )

!*****************************************************************************80
!
!! GAMMA_INC_TEST tests R4_GAMIC and R8_GAMIC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_gamic
  real ( kind = 8 ) r8_gamic
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_INC_TEST:'
  write ( *, '(a)' ) '  Test GAMMA_INC_VALUES, R4_GAMIC, R8_GAMIC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             A               X     GAMIC(A,X)'
  write ( *, '(a)' ) '                                R4_GAMIC(A,X)         Diff'
  write ( *, '(a)' ) '                                R8_GAMIC(A,X)         Diff'

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_gamic ( real ( a ), real ( x, kind = 4 ) )
    fx3 = r8_gamic ( a, x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, x, fx1
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( fx1 - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine gamma_inc_tricomi_test ( )

!*****************************************************************************80
!
!! GAMMA_INC_TRICOMI_TEST tests R4_GAMIT and R8_GAMIT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_gamit
  real ( kind = 8 ) r8_gamit
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_INC_TRICOMI_TEST:'
  write ( *, '(a)' ) '  Test GAMMA_INC_TRICOMI_VALUES, R4_GAMIT, R8_GAMIT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             A               X     GAMIT(A,X)'
  write ( *, '(a)' ) '                                R4_GAMIT(A,X)         Diff'
  write ( *, '(a)' ) '                                R8_GAMIT(A,X)         Diff'

  n_data = 0

  do

    call gamma_inc_tricomi_values ( n_data, a, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_gamit ( real ( a ), real ( x, kind = 4 ) )
    fx3 = r8_gamit ( a, x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, x, fx1
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( fx1 - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine int_test ( )

!*****************************************************************************80
!
!! INT_TEST tests R4_INT and R8_INT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_int
  real ( kind = 8 ) r8_int
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INT_TEST:'
  write ( *, '(a)' ) '  Test INT_VALUES, R4_INT, R8_INT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         INT(X)'
  write ( *, '(a)' ) '                    R4_INT(X)         Diff'
  write ( *, '(a)' ) '                    R8_INT(X)         Diff'

  n_data = 0

  do

    call int_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_int ( real ( x, kind = 4 ) )
    fx3 = r8_int ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine lbeta_test ( )

!*****************************************************************************80
!
!! LBETA_TEST tests R4_LBETA and R8_LBETA.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_lbeta
  real ( kind = 8 ) r8_lbeta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LBETA_TEST:'
  write ( *, '(a)' ) '  Test BETA_LOG_VALUES, R4_LBETA, R8_LBETA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             A               B     LBETA(A,B)'
  write ( *, '(a)' ) '                                R4_LBETA(A,B)         Diff'
  write ( *, '(a)' ) '                                R8_LBETA(A,B)         Diff'

  n_data = 0

  do

    call beta_log_values ( n_data, a, b, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_lbeta ( real ( a ), real ( b ) )
    fx3 = r8_lbeta ( a, b )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, b, fx1
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( fx1 - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine li_test ( )

!*****************************************************************************80
!
!! LI_TEST tests R4_LI and R8_LI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_li
  real ( kind = 8 ) r8_li
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LI_TEST:'
  write ( *, '(a)' ) '  Test LOGARITHMIC_INTEGRAL_VALUES, R4_LI, R8_LI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X           LI(X)'
  write ( *, '(a)' ) '                      R4_LI(X)        Diff'
  write ( *, '(a)' ) '                      R8_LI(X)        Diff'

  n_data = 0

  do

    call logarithmic_integral_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_li ( real ( x, kind = 4 ) )
    fx3 = r8_li ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine lngam_test ( )

!*****************************************************************************80
!
!! LNGAM_TEST tests R4_LNGAM and R8_LNGAM.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_lngam
  real ( kind = 8 ) r8_lngam
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LNGAM_TEST:'
  write ( *, '(a)' ) '  Test GAMMA_LOG_VALUES, R4_LNGAM, R8_LNGAM.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X        LNGAM(X)'
  write ( *, '(a)' ) '                   R4_LNGAM(X)        Diff'
  write ( *, '(a)' ) '                   R8_LNGAM(X)        Diff'

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_lngam ( real ( x, kind = 4 ) )
    fx3 = r8_lngam ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine log_test ( )

!*****************************************************************************80
!
!! LOG_TEST tests R4_LOG and R8_LOG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_log
  real ( kind = 8 ) r8_log
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOG_TEST:'
  write ( *, '(a)' ) '  Test LOG_VALUES, R4_LOG, R8_LOG.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X          LOG(X)'
  write ( *, '(a)' ) '                     R4_LOG(X)        Diff'
  write ( *, '(a)' ) '                     R8_LOG(X)        Diff'

  n_data = 0

  do

    call log_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_log ( real ( x, kind = 4 ) )
    fx3 = r8_log ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine log10_test ( )

!*****************************************************************************80
!
!! LOG10_TEST tests R4_LOG10 and R8_LOG10.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_log10
  real ( kind = 8 ) r8_log10
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOG10_TEST:'
  write ( *, '(a)' ) '  Test LOG10_VALUES, R4_LOG10, R8_LOG10.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X        LOG10(X)'
  write ( *, '(a)' ) '                   R4_LOG10(X)        Diff'
  write ( *, '(a)' ) '                   R8_LOG10(X)        Diff'

  n_data = 0

  do

    call log10_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_log10 ( real ( x, kind = 4 ) )
    fx3 = r8_log10 ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine poch_test ( )

!*****************************************************************************80
!
!! POCH_TEST tests R4_POCH and R8_POCH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_poch
  real ( kind = 8 ) r8_poch
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POCH_TEST:'
  write ( *, '(a)' ) '  Test POCHHAMMER_VALUES, R4_POCH, R8_POCH.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             A               X      POCH(A,X)'
  write ( *, '(a)' ) '                                 R4_POCH(A,X)         Diff'
  write ( *, '(a)' ) '                                 R8_POCH(A,X)         Diff'

  n_data = 0

  do

    call pochhammer_values ( n_data, a, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_poch ( real ( a ), real ( x, kind = 4 ) )
    fx3 = r8_poch ( a, x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.4,2x,g14.6)' )  a, x, fx1
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx2, abs ( fx1 - fx2 )
    write ( *, '(2x,14x,2x,14x,2x,g14.6,2x,g14.6)' ) fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine psi_test ( )

!*****************************************************************************80
!
!! PSI_TEST tests R4_PSI and R8_PSI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_psi
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PSI_TEST:'
  write ( *, '(a)' ) '  Test PSI_VALUES, R4_PSI, R8_PSI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         PSI(X)'
  write ( *, '(a)' ) '                    R4_PSI(X)         Diff'
  write ( *, '(a)' ) '                    R8_PSI(X)         Diff'

  n_data = 0

  do

    call psi_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_psi ( real ( x, kind = 4 ) )
    fx3 = r8_psi ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine rand_test ( )

!*****************************************************************************80
!
!! RAND_TEST tests R4_RAND.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    02 October 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 4 ) average
  integer ( kind = 4 ) i
  integer ( kind = 4 ), save, dimension ( 7 ) :: i_value = (/ &
    1, 2, 3, 4, 10, 100, 1000 /)
  integer ( kind = 4 ) k
  real ( kind = 4 ) r
  real ( kind = 4 ), save, dimension ( 7 ) :: r_value = (/ &
    0.0004127026E+00, &
    0.6750836372E+00, &
    0.1614754200E+00, &
    0.9086198807E+00, &
    0.5527787209E+00, &
    0.3600893021E+00, &
    0.2176990509E+00 /)
  real ( kind = 4 ) r4_rand
  real ( kind = 4 ) variance

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RAND_TEST:'
  write ( *, '(a)' ) '  Test R4_RAND.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             I         R4_RAND        Expected'
  write ( *, '(a)' ) ' '

  k = 1

  do i = 1, 1000

    r = r4_rand ( 0.0E+00 )

    if ( i == i_value(k) ) then
      write ( *, '(2x,i14,2x,g14.6,2x,g14.6)' ) i, r, r_value(k)
      k = k + 1
    end if

  end do

  average = 0.0E+00
  do i = 1, 1000000
    r = r4_rand ( 0.0E+00 )
    average = average + r
  end do
  average = average / 1000000.0E+00
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,2x,g14.6)' ) '       Average =  ', average, 0.5E+00

  variance = 0.0E+00
  do i = 1, 1000000
    r = r4_rand ( 0.0 )
    variance = variance + ( r - average )**2
  end do
  variance = variance / 1000000.0E+00
  write ( *, '(a,g14.6,2x,g14.6)' ) &
     '       Variance = ', variance, 1.0E+00 / 12.0E+00

  return
end
subroutine shi_test ( )

!*****************************************************************************80
!
!! SHI_TEST tests R4_SHI and R8_SHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_shi
  real ( kind = 8 ) r8_shi
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHI_TEST:'
  write ( *, '(a)' ) '  Test SHI_VALUES, R4_SHI, R8_SHI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X          SHI(X)'
  write ( *, '(a)' ) '                     R4_SHI(X)        Diff'
  write ( *, '(a)' ) '                     R8_SHI(X)        Diff'

  n_data = 0

  do

    call shi_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_shi ( real ( x, kind = 4 ) )
    fx3 = r8_shi ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine si_test ( )

!*****************************************************************************80
!
!! SI_TEST tests R4_SI and R8_SI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_si
  real ( kind = 8 ) r8_si
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SI_TEST:'
  write ( *, '(a)' ) '  Test SI_VALUES, R4_SI, R8_SI.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X           SI(X)'
  write ( *, '(a)' ) '                      R4_SI(X)        Diff'
  write ( *, '(a)' ) '                      R8_SI(X)        Diff'

  n_data = 0

  do

    call si_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_si ( real ( x, kind = 4 ) )
    fx3 = r8_si ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine sin_test ( )

!*****************************************************************************80
!
!! SIN_TEST tests R4_SIN and R8_SIN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_sin
  real ( kind = 8 ) r8_sin
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIN_TEST:'
  write ( *, '(a)' ) '  Test SIN_VALUES, R4_SIN, R8_SIN.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         SIN(X)'
  write ( *, '(a)' ) '                    R4_SIN(X)         Diff'
  write ( *, '(a)' ) '                    R8_SIN(X)         Diff'

  n_data = 0

  do

    call sin_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_sin ( real ( x, kind = 4 ) )
    fx3 = r8_sin ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine sin_deg_test ( )

!*****************************************************************************80
!
!! SIN_DEG_TEST tests R4_SIN_DEG and R8_SIN_DEG.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_sin_deg
  real ( kind = 8 ) r8_sin_deg
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIN_DEG_TEST:'
  write ( *, '(a)' ) '  Test SIN_DEGREE_VALUES, R4_SIN_DEG, R8_SIN_DEG.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X     SIN_DEG(X)'
  write ( *, '(a)' ) '                R4_SIN_DEG(X)         Diff'
  write ( *, '(a)' ) '                R8_SIN_DEG(X)         Diff'

  n_data = 0

  do

    call sin_degree_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_sin_deg ( real ( x, kind = 4 ) )
    fx3 = r8_sin_deg ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine sinh_test ( )

!*****************************************************************************80
!
!! SINH_TEST tests R4_SINH and R8_SINH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_sinh
  real ( kind = 8 ) r8_sinh
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINH_TEST:'
  write ( *, '(a)' ) '  Test SINH_VALUES, R4_SINH, R8_SINH.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         SINH(X)'
  write ( *, '(a)' ) '                    R4_SINH(X)        Diff'
  write ( *, '(a)' ) '                    R8_SINH(X)        Diff'

  n_data = 0

  do

    call sinh_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_sinh ( real ( x, kind = 4 ) )
    fx3 = r8_sinh ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine spence_test ( )

!*****************************************************************************80
!
!! SPENCE_TEST tests R4_SPENCE and R8_SPENCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_spence
  real ( kind = 8 ) r8_spence
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPENCE_TEST:'
  write ( *, '(a)' ) '  Test DILOGARITHM_VALUES, R4_SPENCE, R8_SPENCE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X       SPENCE(X)'
  write ( *, '(a)' ) '                  R4_SPENCE(X)        Diff'
  write ( *, '(a)' ) '                  R8_SPENCE(X)        Diff'

  n_data = 0

  do

    call dilogarithm_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_spence ( real ( x, kind = 4 ) )
    fx3 = r8_spence ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine sqrt_test ( )

!*****************************************************************************80
!
!! SQRT_TEST tests R4_SQRT and R8_SQRT.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_sqrt
  real ( kind = 8 ) r8_sqrt
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SQRT_TEST:'
  write ( *, '(a)' ) '  Test SQRT_VALUES, R4_SQRT, R8_SQRT.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         SQRT(X)'
  write ( *, '(a)' ) '                    R4_SQRT(X)        Diff'
  write ( *, '(a)' ) '                    R8_SQRT(X)        Diff'

  n_data = 0

  do

    call sqrt_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_sqrt ( real ( x, kind = 4 ) )
    fx3 = r8_sqrt ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine tan_test ( )

!*****************************************************************************80
!
!! TAN_TEST tests R4_TAN and R8_TAN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_tan
  real ( kind = 8 ) r8_tan
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TAN_TEST:'
  write ( *, '(a)' ) '  Test TAN_VALUES, R4_TAN, R8_TAN.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         TAN(X)'
  write ( *, '(a)' ) '                    R4_TAN(X)         Diff'
  write ( *, '(a)' ) '                    R8_TAN(X)         Diff'

  n_data = 0

  do

    call tan_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_tan ( real ( x, kind = 4 ) )
    fx3 = r8_tan ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
subroutine tanh_test ( )

!*****************************************************************************80
!
!! TANH_TEST tests R4_TANH and R8_TANH.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) diff
  real ( kind = 8 ) fx1
  real ( kind = 4 ) fx2
  real ( kind = 8 ) fx3
  integer ( kind = 4 ) n_data
  real ( kind = 4 ) r4_tanh
  real ( kind = 8 ) r8_tanh
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TANH_TEST:'
  write ( *, '(a)' ) '  Test TANH_VALUES, R4_TANH, R8_TANH.'
  write ( *, '(a)' ) '  TANH_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             X         TANH(X)'
  write ( *, '(a)' ) '                    R4_TANH(X)         Diff'
  write ( *, '(a)' ) '                    R8_TANH(X)         Diff'

  n_data = 0

  do

    call tanh_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_tanh ( real ( x, kind = 4 ) )
    fx3 = r8_tanh ( x )

    write ( *, '(a)' ) ' '
    write ( *, '(2x,f14.4,2x,g14.6)' )     x, fx1
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx2, abs ( fx1 - fx2 )
    write ( *, '(16x,2x,g14.6,2x,g14.6)' )    fx3, abs ( fx1 - fx3 )

  end do

  return
end
