program main

!*****************************************************************************80
!
!! MAIN is the main program for SPECFUN_PRB2.
!
!  Discussion:
!
!    SPECFUN_PRB2 calls sample problems for the SPECFUN library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPECFUN_PRB2'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SPECFUN library.'

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
  write ( *, '(a)' ) 'SPECFUN_PRB2'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 checks BESI0 against BESSEL_I0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besi0
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  BESI0 computes the Bessel I0 function.'
  write ( *, '(a)' ) '  BESSEL_I0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESI0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besi0 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 checks RIBESL against BESSEL_I0_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nb
  parameter ( nb = 1 )

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(nb)
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  RIBESL returns values of Bessel functsions'
  write ( *, '(a)' ) '  of non-integer order.'
  write ( *, '(a)' ) '  BESSEL_I0_SPHERICAL_VALUES returns values'
  write ( *, '(a)' ) '  of the spherical Bessel i0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X                      FX                        FX2'
  write ( *, '(a,a)' ) &
    '                                 (table)', &
    '                   (RIBESL)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    alpha = 0.5D+00
    ize = 1

    call ribesl ( x, alpha, nb, ize, b, ncalc )

    fx2 = sqrt ( 0.5D+00 * pi / x ) * b(1)

    write ( *, '(2x,f24.16,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 checks BESI0 against BESSEL_I0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besi1
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  BESI1 computes the Bessel I1 function.'
  write ( *, '(a)' ) '  BESSEL_I1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESI1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besi1 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 checks RIBESL against BESSEL_I1_SPHERICAL_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nb
  parameter ( nb = 2 )

  real ( kind = 8 ) alpha
  real ( kind = 8 ) b(nb)
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) ize
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  RIBESL returns values of Bessel functsions'
  write ( *, '(a)' ) '  of non-integer order.'
  write ( *, '(a)' ) '  BESSEL_I1_SPHERICAL_VALUES returns values'
  write ( *, '(a)' ) '  of the spherical Bessel i1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X                      FX                        FX2'
  write ( *, '(a,a)' ) &
    '                                 (table)', &
    '                   (RIBESL)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i1_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    alpha = 0.5D+00
    ize = 1

    call ribesl ( x, alpha, nb, ize, b, ncalc )

    fx2 = sqrt ( 0.5D+00 * pi / x ) * b(2)

    write ( *, '(2x,f24.16,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 checks RIBESL against BESSEL_IX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) nb_max
  parameter ( nb_max = 10 )

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
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  RIBESL computes values of Bessel functions'
  write ( *, '(a)' ) '  of NONINTEGER order.'
  write ( *, '(a)' ) '  BESSEL_IX_VALUES returns selected values'
  write ( *, '(a)' ) &
    '  of the Bessel function In for NONINTEGER order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      ALPHA         X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                                  (table)', &
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

    write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      alpha, x, fx, fx2

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 checks BESJ0 against BESSEL_J0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besj0
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  BESJ0 computes the Bessel J0 function.'
  write ( *, '(a)' ) '  BESSEL_J0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESJ0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besj0 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 checks BESJ1 against BESSEL_J1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besj1
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  BESJ1 computes the Bessel J1 function.'
  write ( *, '(a)' ) '  BESSEL_J1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESJ1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besj1 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 checks RJBESL against BESSEL_JX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
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
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  RJBESL computes values of Bessel functions'
  write ( *, '(a)' ) '  of NONINTEGER order.'
  write ( *, '(a)' ) '  BESSEL_JX_VALUES returns selected values'
  write ( *, '(a)' ) &
    '  of the Bessel function Jn for NONINTEGER order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      ALPHA         X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                                  (table)', &
    '                   (RJBESL)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jx_values ( n_data, alpha, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    nb = int ( alpha ) + 1
    if ( nb_max < nb ) then
      write ( *, * ) '  [Skipping calculation, NB_MAX too small.]'
      cycle
    end if

    alpha_frac = alpha - real ( int ( alpha ), kind = 8 )

    call rjbesl ( x, alpha_frac, nb, b, ncalc )

    fx2 = b(nb)

    write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      alpha, x, fx, fx2

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 checks BESK0 against BESSEL_K0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besk0
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  BESK0 computes the Bessel K0 function.'
  write ( *, '(a)' ) '  BESSEL_K0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESK0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besk0 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 checks BESK1 against BESSEL_K1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besk1
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10:'
  write ( *, '(a)' ) '  BESK1 computes the Bessel K1 function.'
  write ( *, '(a)' ) '  BESSEL_K1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESK1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besk1 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 checks RKBESL against BESSEL_KX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
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
  write ( *, '(a)' ) 'TEST11:'
  write ( *, '(a)' ) '  RKBESL computes values of Bessel functions'
  write ( *, '(a)' ) '  of NONINTEGER order.'
  write ( *, '(a)' ) '  BESSEL_KX_VALUES returns selected values'
  write ( *, '(a)' ) &
    '  of the Bessel function Kn for NONINTEGER order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      ALPHA         X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                                  (table)', &
    '                   (RKBESL)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_kx_values ( n_data, alpha, x, fx )

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

    call rkbesl ( x, alpha_frac, nb, ize, b, ncalc )

    fx2 = b(nb)

    write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      alpha, x, fx, fx2

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 checks BESY0 against BESSEL_Y0_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besy0
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12:'
  write ( *, '(a)' ) '  BESY0 computes the Bessel Y0 function.'
  write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESY0)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besy0 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 checks BESY1 against BESSEL_Y1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) besy1
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13:'
  write ( *, '(a)' ) '  BESY1 computes the Bessel Y1 function.'
  write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (BESY1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = besy1 ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 checks RYBESL against BESSEL_YX_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
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
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) ncalc
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14:'
  write ( *, '(a)' ) '  RYBESL computes values of Bessel functions'
  write ( *, '(a)' ) '  of NONINTEGER order.'
  write ( *, '(a)' ) '  BESSEL_YX_VALUES returns selected values'
  write ( *, '(a)' ) &
    '  of the Bessel function Yn for NONINTEGER order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      ALPHA         X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                                  (table)', &
    '                   (RYBESL)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_yx_values ( n_data, alpha, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    nb = int ( alpha ) + 1
    if ( nb_max < nb ) then
      write ( *, * ) '  [Skipping calculation, NB_MAX too small.]'
      cycle
    end if

    alpha_frac = alpha - real ( int ( alpha ), kind = 8 )

    call rybesl ( x, alpha_frac, nb, b, ncalc )

    fx2 = b(nb)

    write ( *, '(2x,f12.8,2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      alpha, x, fx, fx2

  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 checks DAW against DAWSON_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) daw
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15:'
  write ( *, '(a)' ) '  DAW computes the Dawson function.'
  write ( *, '(a)' ) '  DAWSON_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DAW)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dawson_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = daw ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 checks EONE against E1_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) eone
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16:'
  write ( *, '(a)' ) '  EONE computes the exponential integral E1.'
  write ( *, '(a)' ) '  E1_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (E1)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call e1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = eone ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 checks EI against EI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ei
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17:'
  write ( *, '(a)' ) '  EI computes the exponential integral Ei.'
  write ( *, '(a)' ) '  EI_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (Ei)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ei_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = ei ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 checks R8_ERF against ERF_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_erf
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18:'
  write ( *, '(a)' ) '  R8_ERF computes the error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DERF)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_erf ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 checks R8_GAMMA against GAMMA_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_gamma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19:'
  write ( *, '(a)' ) &
    '  DGAMMA computes the gamma function.'
  write ( *, '(a)' ) '  GAMMA_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2                    DIFF'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (R8_GAMMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_gamma ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 checks DLGAMA against GAMMA_LOG_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) dlgama
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20:'
  write ( *, '(a)' ) &
    '  DLGAMA computes the log of the gamma function.'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             FX', &
    '                        FX2'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (DLGAMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dlgama ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16)' ) &
      x, fx, fx2

  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests MACHAR.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) eps
  real ( kind = 8 ) epsneg
  integer ( kind = 4 ) ibeta
  integer ( kind = 4 ) iexp
  integer ( kind = 4 ) irnd
  integer ( kind = 4 ) it
  integer ( kind = 4 ) machep
  integer ( kind = 4 ) maxexp
  integer ( kind = 4 ) minexp
  integer ( kind = 4 ) negep
  integer ( kind = 4 ) ngrd
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  MACHAR computes double'
  write ( *, '(a)' ) '  precision machine constants.'

  call machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, &
    minexp, maxexp, eps, epsneg, xmin, xmax )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  IBETA is the internal base for machine arithmetic.'
  write ( *, '(a,i8)' ) '    IBETA =  ', ibeta
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  IT is the number of digits, base IBETA, in the'
  write ( *, '(a)' ) '  floating point significand.'
  write ( *, '(a,i8)' ) '    IT =     ', it
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  IRND reports on floating point addition rounding:'
  write ( *, '(a)' ) '  0, for chopping;'
  write ( *, '(a)' ) '  1, for non-IEEE rounding;'
  write ( *, '(a)' ) '  2, for IEEE rounding;'
  write ( *, '(a)' ) '  3, for chopping with partial underflow;'
  write ( *, '(a)' ) &
    '  4, for non-IEEE rounding with partial underflow.'
  write ( *, '(a)' ) &
    '  5, for IEEE rounding with partial underflow.'
  write ( *, '(a,i8)' ) '    IRND =   ', irnd
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  NGRD is the number of guard digits for floating point'
  write ( *, '(a)' ) '  multiplication with truncating arithmetic.'
  write ( *, '(a,i8)' ) '    NGRD =   ', ngrd
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  MACHEP is the largest negative integer such that'
  write ( *, '(a)' ) '  1.0 < 1.0 + BETA**MACHEP.'
  write ( *, '(a,i8)' ) '    MACHEP = ', machep
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  NEGEPS is the largest negative integer such that'
  write ( *, '(a)' ) '  1.0 - BETA**NEGEPS < 1.0:'
  write ( *, '(a,i8)' ) '    NEGEP =  ', negep
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  IEXP is the number of bits reserved for the exponent'
  write ( *, '(a)' ) '  of a floating point number:'
  write ( *, '(a,i8)' ) '    IEXP =   ', iexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  MINEXP is the most negative power of BETA such that'
  write ( *, '(a)' ) '  BETA**MINEXP is positive and normalized.'
  write ( *, '(a,i8)' ) '    MINEXP = ', minexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  MAXEXP is the smallest positive power of BETA that'
  write ( *, '(a)' ) '  overflows:'
  write ( *, '(a,i8)' ) '    MAXEXP = ', maxexp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  EPS is a small positive floating point number'
  write ( *, '(a)' ) '  such that 1.0 < 1.0 + EPS.'
  write ( *, '(a,e25.13)' ) '    EPS    = ', eps
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  EPSNEG is a small positive floating point number'
  write ( *, '(a)' ) '  such that 1.0 - EPSNEG < 1.0.'
  write ( *, '(a,e25.13)' ) '    EPSNEG = ', epsneg
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  XMIN is the smallest positive normalized floating'
  write ( *, '(a)' ) '  point power of the radix:'
  write ( *, '(a,e25.13)' ) '    XMIN =   ', xmin
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  XMAX is the largest finite floating point number:'
  write ( *, '(a,e25.13)' ) '    XMAX   = ', xmax

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 checks R8_PSI against PSI_VALUES.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) r8_psi
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22:'
  write ( *, '(a)' ) '  R8_PSI computes the PSI function.'
  write ( *, '(a)' ) '  PSI_VALUES returns selected values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      X             PSI', &
    '                       PSI                    DIFF'
  write ( *, '(a,a)' ) &
    '                   (table)', &
    '                   (R8_PSI)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_psi ( x )

    write ( *, '(2x,f12.8,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test23 ( )

!*********************************************************************80
!
!! TEST23 tests REN.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    04 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: sample_num = 1000

  integer ( kind = 4 ) i
  real ( kind = 8 ) mean
  real ( kind = 8 ) ren
  integer ( kind = 4 ) seed
  real ( kind = 8 ) variance
  real ( kind = 8 ) x(sample_num)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  REN is a random number generator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  It cannot be controlled by an external'
  write ( *, '(a)' ) '  seed value, although it has an argument'
  write ( *, '(a)' ) '  that "looks" like a seed.'

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  10 sample values:'
  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(g14.6)' ) ren ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Expected mean =     ', 0.5D+00
  write ( *, '(a,g14.6)' ) '  Expected variance = ', &
    1.0D+00 / 12.0D+00

  mean = 0.0D+00
  xmax = -1000.0D+00
  xmin = +1000.0D+00

  do i = 1, sample_num
    x(i) = ren ( seed )
    mean = mean + x(i)
    xmax = max ( xmax, x(i) )
    xmin = min ( xmin, x(i) )
  end do
  mean = mean / real ( sample_num, kind = 8 )

  variance = 0.0D+00
  do i = 1, sample_num
    variance = variance + ( x(i) - mean )**2
  end do
  variance = variance / real ( sample_num - 1, kind = 8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'    ) '  Sample size =     ', sample_num
  write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
  write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
  write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
  write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

  return
end
subroutine bessel_i0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I0_VALUES returns some values of the I0 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    The modified Bessel function I0(Z) corresponds to N = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.1000000000000000D+01, &
    0.1010025027795146D+01, &
    0.1040401782229341D+01, &
    0.1092045364317340D+01, &
    0.1166514922869803D+01, &
    0.1266065877752008D+01, &
    0.1393725584134064D+01, &
    0.1553395099731217D+01, &
    0.1749980639738909D+01, &
    0.1989559356618051D+01, &
    0.2279585302336067D+01, &
    0.3289839144050123D+01, &
    0.4880792585865024D+01, &
    0.7378203432225480D+01, &
    0.1130192195213633D+02, &
    0.1748117185560928D+02, &
    0.2723987182360445D+02, &
    0.6723440697647798D+02, &
    0.4275641157218048D+03, &
    0.2815716628466254D+04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    0.10D+01, &
    0.12D+01, &
    0.14D+01, &
    0.16D+01, &
    0.18D+01, &
    0.20D+01, &
    0.25D+01, &
    0.30D+01, &
    0.35D+01, &
    0.40D+01, &
    0.45D+01, &
    0.50D+01, &
    0.60D+01, &
    0.80D+01, &
    0.10D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_i0_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I0_SPHERICAL_VALUES returns some values of the Spherical Bessel function i0.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselI[1/2,x]
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    1.001667500198440D+00, &
    1.006680012705470D+00, &
    1.026880814507039D+00, &
    1.061089303580402D+00, &
    1.110132477734529D+00, &
    1.175201193643801D+00, &
    1.257884462843477D+00, &
    1.360215358179667D+00, &
    1.484729970750144D+00, &
    1.634541271164267D+00, &
    1.813430203923509D+00, &
    2.025956895698133D+00, &
    2.277595505698373D+00, &
    2.574897010920645D+00, &
    2.925685126512827D+00, &
    3.339291642469967D+00, &
    3.826838748926716D+00, &
    4.401577467270101D+00, &
    5.079293155726485D+00, &
    5.878791279137455D+00, &
    6.822479299281938D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_i1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I1_VALUES returns some values of the I1 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1005008340281251D+00, &
    0.2040267557335706D+00, &
    0.3137040256049221D+00, &
    0.4328648026206398D+00, &
    0.5651591039924850D+00, &
    0.7146779415526431D+00, &
    0.8860919814143274D+00, &
    0.1084810635129880D+01, &
    0.1317167230391899D+01, &
    0.1590636854637329D+01, &
    0.2516716245288698D+01, &
    0.3953370217402609D+01, &
    0.6205834922258365D+01, &
    0.9759465153704450D+01, &
    0.1538922275373592D+02, &
    0.2433564214245053D+02, &
    0.6134193677764024D+02, &
    0.3998731367825601D+03, &
    0.2670988303701255D+04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.00D+00, &
    0.20D+00, &
    0.40D+00, &
    0.60D+00, &
    0.80D+00, &
    0.10D+01, &
    0.12D+01, &
    0.14D+01, &
    0.16D+01, &
    0.18D+01, &
    0.20D+01, &
    0.25D+01, &
    0.30D+01, &
    0.35D+01, &
    0.40D+01, &
    0.45D+01, &
    0.50D+01, &
    0.60D+01, &
    0.80D+01, &
    0.10D+02 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_i1_spherical_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_I1_SPHERICAL_VALUES returns some values of the Spherical Bessel function i1.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi/(2*x)] * BesselI[3/2,x]
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.03336667857363341D+00, &
    0.06693371456802954D+00, &
    0.1354788933285401D+00, &
    0.2072931911031093D+00, &
    0.2841280857128948D+00, &
    0.3678794411714423D+00, &
    0.4606425870674146D+00, &
    0.5647736480096238D+00, &
    0.6829590627779635D+00, &
    0.8182955028627777D+00, &
    0.9743827435800610D+00, &
    1.155432469636406D+00, &
    1.366396525527973D+00, &
    1.613118767572064D+00, &
    1.902515460838681D+00, &
    2.242790117769266D+00, &
    2.643689828630357D+00, &
    3.116811526884873D+00, &
    3.675968313148932D+00, &
    4.337627987747642D+00, &
    5.121438384183637D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.1D+00, &
    0.2D+00, &
    0.4D+00, &
    0.6D+00, &
    0.8D+00, &
    1.0D+00, &
    1.2D+00, &
    1.4D+00, &
    1.6D+00, &
    1.8D+00, &
    2.0D+00, &
    2.2D+00, &
    2.4D+00, &
    2.6D+00, &
    2.8D+00, &
    3.0D+00, &
    3.2D+00, &
    3.4D+00, &
    3.6D+00, &
    3.8D+00, &
    4.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_ix_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_IX_VALUES returns some values of the Ix Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function In is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "In" by "Ix".
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselI[n,x]
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.3592084175833614D+00,  &
    0.9376748882454876D+00,  &
    2.046236863089055D+00,   &
    3.053093538196718D+00,   &
    4.614822903407601D+00,   &
    26.47754749755907D+00,   &
    2778.784603874571D+00,   &
    4.327974627242893D+07,   &
    0.2935253263474798D+00,  &
    1.099473188633110D+00,   &
    21.18444226479414D+00,   &
    2500.906154942118D+00,   &
    2.866653715931464D+20,   &
    0.05709890920304825D+00, &
    0.3970270801393905D+00,  &
    13.76688213868258D+00,   &
    2028.512757391936D+00,   &
    2.753157630035402D+20,   &
    0.4139416015642352D+00,  &
    1.340196758982897D+00,   &
    22.85715510364670D+00,   &
    2593.006763432002D+00,   &
    2.886630075077766D+20,   &
    0.03590910483251082D+00, &
    0.2931108636266483D+00,  &
    11.99397010023068D+00,   &
    1894.575731562383D+00,   &
    2.716911375760483D+20 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_j0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J0_VALUES returns some values of the J0 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1775967713143383D+00, &
    -0.3971498098638474D+00, &
    -0.2600519549019334D+00, &
     0.2238907791412357D+00, &
     0.7651976865579666D+00, &
     0.1000000000000000D+01, &
     0.7651976865579666D+00, &
     0.2238907791412357D+00, &
    -0.2600519549019334D+00, &
    -0.3971498098638474D+00, &
    -0.1775967713143383D+00, &
     0.1506452572509969D+00, &
     0.3000792705195556D+00, &
     0.1716508071375539D+00, &
    -0.9033361118287613D-01, &
    -0.2459357644513483D+00, &
    -0.1711903004071961D+00, &
     0.4768931079683354D-01, &
     0.2069261023770678D+00, &
     0.1710734761104587D+00, &
    -0.1422447282678077D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -5.0D+00, &
    -4.0D+00, &
    -3.0D+00, &
    -2.0D+00, &
    -1.0D+00, &
     0.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_j1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_J1_VALUES returns some values of the J1 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.3275791375914652D+00, &
     0.6604332802354914D-01, &
    -0.3390589585259365D+00, &
    -0.5767248077568734D+00, &
    -0.4400505857449335D+00, &
     0.0000000000000000D+00, &
     0.4400505857449335D+00, &
     0.5767248077568734D+00, &
     0.3390589585259365D+00, &
    -0.6604332802354914D-01, &
    -0.3275791375914652D+00, &
    -0.2766838581275656D+00, &
    -0.4682823482345833D-02, &
     0.2346363468539146D+00, &
     0.2453117865733253D+00, &
     0.4347274616886144D-01, &
    -0.1767852989567215D+00, &
    -0.2234471044906276D+00, &
    -0.7031805212177837D-01, &
     0.1333751546987933D+00, &
     0.2051040386135228D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -5.0D+00, &
    -4.0D+00, &
    -3.0D+00, &
    -2.0D+00, &
    -1.0D+00, &
     0.0D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)
  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_jx_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_JX_VALUES returns some values of the Jx Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function Jn is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "Jn" by "Jx".
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselJ[n,x]
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
        0.3544507442114011D+00, &
        0.6713967071418031D+00, &
        0.5130161365618278D+00, &
        0.3020049060623657D+00, &
        0.06500818287737578D+00, &
       -0.3421679847981618D+00, &
       -0.1372637357550505D+00, &
        0.1628807638550299D+00, &
        0.2402978391234270D+00, &
        0.4912937786871623D+00, &
       -0.1696513061447408D+00, &
        0.1979824927558931D+00, &
       -0.1094768729883180D+00, &
        0.04949681022847794D+00, &
        0.2239245314689158D+00, &
        0.2403772011113174D+00, &
        0.1966584835818184D+00, &
        0.02303721950962553D+00, &
        0.3314145508558904D+00, &
        0.5461734240402840D+00, &
       -0.2616584152094124D+00, &
        0.1296035513791289D+00, &
       -0.1117432171933552D+00, &
        0.03142623570527935D+00, &
        0.1717922192746527D+00, &
        0.3126634069544786D+00, &
        0.1340289119304364D+00, &
        0.06235967135106445D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_k0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_K0_VALUES returns some values of the K0 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    The modified Bessel function K0(Z) corresponds to N = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselK[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.2427069024702017D+01, &
    0.1752703855528146D+01, &
    0.1114529134524434D+01, &
    0.7775220919047293D+00, &
    0.5653471052658957D+00, &
    0.4210244382407083D+00, &
    0.3185082202865936D+00, &
    0.2436550611815419D+00, &
    0.1879547519693323D+00, &
    0.1459314004898280D+00, &
    0.1138938727495334D+00, &
    0.6234755320036619D-01, &
    0.3473950438627925D-01, &
    0.1959889717036849D-01, &
    0.1115967608585302D-01, &
    0.6399857243233975D-02, &
    0.3691098334042594D-02, &
    0.1243994328013123D-02, &
    0.1464707052228154D-03, &
    0.1778006231616765D-04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     0.2D+00, &
     0.4D+00, &
     0.6D+00, &
     0.8D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00, &
     5.0D+00, &
     6.0D+00, &
     8.0D+00, &
    10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_k1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_K1_VALUES returns some values of the K1 Bessel function.
!
!  Discussion:
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    The modified Bessel function K1(Z) corresponds to N = 1.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselK[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.9853844780870606D+01, &
    0.4775972543220472D+01, &
    0.2184354424732687D+01, &
    0.1302834939763502D+01, &
    0.8617816344721803D+00, &
    0.6019072301972346D+00, &
    0.4345923910607150D+00, &
    0.3208359022298758D+00, &
    0.2406339113576119D+00, &
    0.1826230998017470D+00, &
    0.1398658818165224D+00, &
    0.7389081634774706D-01, &
    0.4015643112819418D-01, &
    0.2223939292592383D-01, &
    0.1248349888726843D-01, &
    0.7078094908968090D-02, &
    0.4044613445452164D-02, &
    0.1343919717735509D-02, &
    0.1553692118050011D-03, &
    0.1864877345382558D-04 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     0.2D+00, &
     0.4D+00, &
     0.6D+00, &
     0.8D+00, &
     1.0D+00, &
     1.2D+00, &
     1.4D+00, &
     1.6D+00, &
     1.8D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     3.5D+00, &
     4.0D+00, &
     4.5D+00, &
     5.0D+00, &
     6.0D+00, &
     8.0D+00, &
    10.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_kx_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_KX_VALUES returns some values of the Kx Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function Kn is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "Kn" by "Kx".
!
!    The modified Bessel functions In(Z) and Kn(Z) are solutions of
!    the differential equation
!
!      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselK[n,x]
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
       2.294489339798475D+00, &
       0.4610685044478946D+00, &
       0.1199377719680614D+00, &
       0.06506594315400999D+00, &
       0.03602598513176459D+00, &
       0.003776613374642883D+00, &
       0.00001799347809370518D+00, &
       5.776373974707445D-10, &
       0.9221370088957891D+00, &
       0.1799066579520922D+00, &
       0.004531936049571459D+00, &
       0.00001979282590307570D+00, &
       3.486992497366216D-23, &
       3.227479531135262D+00, &
       0.3897977588961997D+00, &
       0.006495775004385758D+00, &
       0.00002393132586462789D+00, &
       3.627839645299048D-23, &
       0.7311451879202114D+00, &
       0.1567475478393932D+00, &
       0.004257389528177461D+00, &
       0.00001915541065869563D+00, &
       3.463337593569306D-23, &
       4.731184839919541D+00, &
       0.4976876225514758D+00, &
       0.007300864610941163D+00, &
       0.00002546421294106458D+00, &
       3.675275677913656D-23 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y0_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y0_VALUES returns some values of the Y0 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[0,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.1534238651350367D+01, &
     0.8825696421567696D-01, &
     0.5103756726497451D+00, &
     0.3768500100127904D+00, &
    -0.1694073932506499D-01, &
    -0.3085176252490338D+00, &
    -0.2881946839815792D+00, &
    -0.2594974396720926D-01, &
     0.2235214893875662D+00, &
     0.2499366982850247D+00, &
     0.5567116728359939D-01, &
    -0.1688473238920795D+00, &
    -0.2252373126343614D+00, &
    -0.7820786452787591D-01, &
     0.1271925685821837D+00, &
     0.2054642960389183D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_y1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! BESSEL_Y1_VALUES returns some values of the Y1 Bessel function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    16 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.6458951094702027D+01, &
    -0.7812128213002887D+00, &
    -0.1070324315409375D+00, &
     0.3246744247918000D+00, &
     0.3979257105571000D+00, &
     0.1478631433912268D+00, &
    -0.1750103443003983D+00, &
    -0.3026672370241849D+00, &
    -0.1580604617312475D+00, &
     0.1043145751967159D+00, &
     0.2490154242069539D+00, &
     0.1637055374149429D+00, &
    -0.5709921826089652D-01, &
    -0.2100814084206935D+00, &
    -0.1666448418561723D+00, &
     0.2107362803687351D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.1D+00, &
     1.0D+00, &
     2.0D+00, &
     3.0D+00, &
     4.0D+00, &
     5.0D+00, &
     6.0D+00, &
     7.0D+00, &
     8.0D+00, &
     9.0D+00, &
    10.0D+00, &
    11.0D+00, &
    12.0D+00, &
    13.0D+00, &
    14.0D+00, &
    15.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine bessel_yx_values ( n_data, nu, x, fx )

!*****************************************************************************80
!
!! BESSEL_YX_VALUES returns some values of the Yx Bessel function.
!
!  Discussion:
!
!    This set of data considers the less common case in which the
!    index of the Bessel function Yn is actually not an integer.
!    We may suggest this case by occasionally replacing the symbol
!    "Yn" by "Yx".
!
!    In Mathematica, the function can be evaluated by:
!
!      BesselY[n,x]
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) NU, the order of the function.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 28

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
       -1.748560416961876D+00, &
       -0.4310988680183761D+00, &
        0.2347857104062485D+00, &
        0.4042783022390569D+00, &
        0.4560488207946332D+00, &
       -0.1012177091851084D+00, &
        0.2117088663313982D+00, &
       -0.07280690478506185D+00, &
       -1.102495575160179D+00, &
       -0.3956232813587035D+00, &
        0.3219244429611401D+00, &
        0.1584346223881903D+00, &
        0.02742813676191382D+00, &
       -2.876387857462161D+00, &
       -0.8282206324443037D+00, &
        0.2943723749617925D+00, &
       -0.1641784796149411D+00, &
        0.1105304445562544D+00, &
       -0.9319659251969881D+00, &
       -0.2609445010948933D+00, &
        0.2492796362185881D+00, &
        0.2174410301416733D+00, &
       -0.01578576650557229D+00, &
       -4.023453301501028D+00, &
       -0.9588998694752389D+00, &
        0.2264260361047367D+00, &
       -0.2193617736566760D+00, &
        0.09413988344515077D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) nu
  real ( kind = 8 ), save, dimension ( n_max ) :: nu_vec = (/ &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    0.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    1.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    2.50D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    1.25D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00, &
    2.75D+00 /)
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.2D+00, &
     1.0D+00, &
     2.0D+00, &
     2.5D+00, &
     3.0D+00, &
     5.0D+00, &
    10.0D+00, &
    20.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00, &
     1.0D+00, &
     2.0D+00, &
     5.0D+00, &
    10.0D+00, &
    50.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    nu = 0.0D+00
    x = 0.0D+00
    fx = 0.0D+00
  else
    nu = nu_vec(n_data)
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine dawson_values ( n_data, x, fx )

!*****************************************************************************80
!
!! DAWSON_VALUES returns some values of Dawson's integral.
!
!  Discussion:
!
!    The definition of Dawson's integral is
!
!      D(X) = exp ( -X * X ) * Integral ( 0 <= Y <= X ) exp ( Y * Y ) dY
!
!    Dawson's integral has a maximum at roughly
!
!      X = 0.9241388730
!
!    In Mathematica, the function can be evaluated by:
!
!      Sqrt[Pi] * Exp[-x^2] * I * Erf[I*x] / 2
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Eric Weisstein,
!    CRC Concise Encyclopedia of Mathematics,
!    CRC Press, 1998.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.9933599239785286D-01, &
    0.1947510333680280D+00, &
    0.2826316650213119D+00, &
    0.3599434819348881D+00, &
    0.4244363835020223D+00, &
    0.4747632036629779D+00, &
    0.5105040575592318D+00, &
    0.5321017070563654D+00, &
    0.5407243187262987D+00, &
    0.5380795069127684D+00, &
    0.5262066799705525D+00, &
    0.5072734964077396D+00, &
    0.4833975173848241D+00, &
    0.4565072375268973D+00, &
    0.4282490710853986D+00, &
    0.3999398943230814D+00, &
    0.3725593489740788D+00, &
    0.3467727691148722D+00, &
    0.3229743193228178D+00, &
    0.3013403889237920D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine e1_values ( n_data, x, fx )

!*****************************************************************************80
!
!! E1_VALUES returns some values of the exponential integral function E1(X).
!
!  Discussion:
!
!    The exponential integral E1(X) is defined by the formula:
!
!      E1(X) = integral ( 1 <= T <= Infinity ) exp ( -X*T ) / T dT
!
!    In Mathematica, the function can be evaluated by:
!
!      ExpIntegralE[1,x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.5597735947761608D+00, &
    0.4543795031894021D+00, &
    0.3737688432335091D+00, &
    0.3105965785455430D+00, &
    0.2601839393259996D+00, &
    0.2193839343955203D+00, &
    0.1859909045360402D+00, &
    0.1584084368514626D+00, &
    0.1354509578491291D+00, &
    0.1162193125713579D+00, &
    0.1000195824066327D+00, &
    0.8630833369753979D-01, &
    0.7465464440125305D-01, &
    0.6471312936386886D-01, &
    0.5620437817453485D-01, &
    0.4890051070806112D-01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine ei_values ( n_data, x, fx )

!*****************************************************************************80
!
!! EI_VALUES returns some values of the exponential integral function EI(X).
!
!  Discussion:
!
!    The exponential integral EI(X) has the formula:
!
!      EI(X) = - integral ( -X <= T <= Infinity ) exp ( -T ) / T dT
!
!    In Mathematica, the function can be evaluated by:
!
!      ExpIntegralEi[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    19 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 16

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.4542199048631736D+00, &
    0.7698812899373594D+00, &
    0.1064907194624291D+01, &
    0.1347396548212326D+01, &
    0.1622811713696867D+01, &
    0.1895117816355937D+01, &
    0.2167378279563403D+01, &
    0.2442092285192652D+01, &
    0.2721398880232024D+01, &
    0.3007207464150646D+01, &
    0.3301285449129798D+01, &
    0.3605319949019469D+01, &
    0.3920963201354904D+01, &
    0.4249867557487934D+01, &
    0.4593713686953585D+01, &
    0.4954234356001890D+01 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine erf_values ( n_data, x, fx )

!*****************************************************************************80
!
!! ERF_VALUES returns some values of the ERF or "error" function.
!
!  Discussion:
!
!    The error function is defined by:
!
!      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
!
!    In Mathematica, the function can be evaluated by:
!
!      Erf[x]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 21

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    0.0000000000000000D+00, &
    0.1124629160182849D+00, &
    0.2227025892104785D+00, &
    0.3286267594591274D+00, &
    0.4283923550466685D+00, &
    0.5204998778130465D+00, &
    0.6038560908479259D+00, &
    0.6778011938374185D+00, &
    0.7421009647076605D+00, &
    0.7969082124228321D+00, &
    0.8427007929497149D+00, &
    0.8802050695740817D+00, &
    0.9103139782296354D+00, &
    0.9340079449406524D+00, &
    0.9522851197626488D+00, &
    0.9661051464753107D+00, &
    0.9763483833446440D+00, &
    0.9837904585907746D+00, &
    0.9890905016357307D+00, &
    0.9927904292352575D+00, &
    0.9953222650189527D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    0.0D+00, &
    0.1D+00, &
    0.2D+00, &
    0.3D+00, &
    0.4D+00, &
    0.5D+00, &
    0.6D+00, &
    0.7D+00, &
    0.8D+00, &
    0.9D+00, &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_VALUES returns some values of the Gamma function.
!
!  Discussion:
!
!    The Gamma function is defined as:
!
!      Gamma(Z) = Integral ( 0 <= T < Infinity) T**(Z-1) exp(-T) dT
!
!    It satisfies the recursion:
!
!      Gamma(X+1) = X * Gamma(X)
!
!    Gamma is undefined for nonpositive integral X.
!    Gamma(0.5) = sqrt(PI)
!    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
!
!    In Mathematica, the function can be evaluated by:
!
!      Gamma[x]
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
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 25

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.3544907701811032D+01, &
    -0.1005871979644108D+03, &
     0.9943258511915060D+02, &
     0.9513507698668732D+01, &
     0.4590843711998803D+01, &
     0.2218159543757688D+01, &
     0.1772453850905516D+01, &
     0.1489192248812817D+01, &
     0.1164229713725303D+01, &
     0.1000000000000000D+01, &
     0.9513507698668732D+00, &
     0.9181687423997606D+00, &
     0.8974706963062772D+00, &
     0.8872638175030753D+00, &
     0.8862269254527580D+00, &
     0.8935153492876903D+00, &
     0.9086387328532904D+00, &
     0.9313837709802427D+00, &
     0.9617658319073874D+00, &
     0.1000000000000000D+01, &
     0.2000000000000000D+01, &
     0.6000000000000000D+01, &
     0.3628800000000000D+06, &
     0.1216451004088320D+18, &
     0.8841761993739702D+31 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    -0.50D+00, &
    -0.01D+00, &
     0.01D+00, &
     0.10D+00, &
     0.20D+00, &
     0.40D+00, &
     0.50D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine gamma_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! GAMMA_LOG_VALUES returns some values of the Log Gamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      Log[Gamma[x]]
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 20

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
     0.1524063822430784D+01, &
     0.7966778177017837D+00, &
     0.3982338580692348D+00, &
     0.1520596783998375D+00, &
     0.0000000000000000D+00, &
    -0.4987244125983972D-01, &
    -0.8537409000331584D-01, &
    -0.1081748095078604D+00, &
    -0.1196129141723712D+00, &
    -0.1207822376352452D+00, &
    -0.1125917656967557D+00, &
    -0.9580769740706586D-01, &
    -0.7108387291437216D-01, &
    -0.3898427592308333D-01, &
    0.00000000000000000D+00, &
    0.69314718055994530D+00, &
    0.17917594692280550D+01, &
    0.12801827480081469D+02, &
    0.39339884187199494D+02, &
    0.71257038967168009D+02 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
     0.20D+00, &
     0.40D+00, &
     0.60D+00, &
     0.80D+00, &
     1.00D+00, &
     1.10D+00, &
     1.20D+00, &
     1.30D+00, &
     1.40D+00, &
     1.50D+00, &
     1.60D+00, &
     1.70D+00, &
     1.80D+00, &
     1.90D+00, &
     2.00D+00, &
     3.00D+00, &
     4.00D+00, &
    10.00D+00, &
    20.00D+00, &
    30.00D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
subroutine psi_values ( n_data, x, fx )

!*****************************************************************************80
!
!! PSI_VALUES returns some values of the Psi or Digamma function.
!
!  Discussion:
!
!    In Mathematica, the function can be evaluated by:
!
!      PolyGamma[x]
!
!    or
!
!      PolyGamma[0,x]
!
!    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
!
!    PSI(1) = -Euler's constant.
!
!    PSI(X+1) = PSI(X) + 1 / X.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    17 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
!    first call.  On each call, the routine increments N_DATA by 1, and
!    returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, real ( kind = 8 ) X, the argument of the function.
!
!    Output, real ( kind = 8 ) FX, the value of the function.
!
  implicit none

  integer ( kind = 4 ), parameter :: n_max = 11

  real ( kind = 8 ) fx
  real ( kind = 8 ), save, dimension ( n_max ) :: fx_vec = (/ &
    -0.5772156649015329D+00, &
    -0.4237549404110768D+00, &
    -0.2890398965921883D+00, &
    -0.1691908888667997D+00, &
    -0.6138454458511615D-01, &
     0.3648997397857652D-01, &
     0.1260474527734763D+00, &
     0.2085478748734940D+00, &
     0.2849914332938615D+00, &
     0.3561841611640597D+00, &
     0.4227843350984671D+00 /)
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ), save, dimension ( n_max ) :: x_vec = (/ &
    1.0D+00, &
    1.1D+00, &
    1.2D+00, &
    1.3D+00, &
    1.4D+00, &
    1.5D+00, &
    1.6D+00, &
    1.7D+00, &
    1.8D+00, &
    1.9D+00, &
    2.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    fx = 0.0D+00
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
