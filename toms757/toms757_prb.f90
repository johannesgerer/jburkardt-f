program main

!*****************************************************************************80
!
!! MAIN is the main program for TOMS757_PRB.
!
!  Discussion:
!
!    TOMS757_PRB tests the routines in TOMS757.
!
!    This program tests the 37 functions in the MISCFUN package.
!    It is a fairly simple code with each function being tested
!    at 20 different arguments.  The code compares the value
!    from the function with a pre-computed value, and produces
!    the absolute and relative errors.
!
!  Author:
!
!    Allan McLeod,
!    Department of Mathematics and Statistics,
!    Paisley University, High Street, Paisley, Scotland, PA12BE
!    macl_ms0@paisley.ac.uk
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS757_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TOMS757 library.'

  call test01
  call test02
  call test03
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09

  call test10
  call test11
  call test12
  call test13
  call test14
  call test15
  call test16
  call test17
  call test18
  call test19

  call test20
  call test21
  call test22
  call test23
  call test24
  call test25
  call test26
  call test27
  call test28
  call test29

  call test30
  call test31
  call test32
  call test33
  call test34
  call test35
  call test36
  call test37
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS757_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01

!*****************************************************************************80
!
!! TEST01 tests ABRAM0.
!
  implicit none

  real ( kind = 8 ) abram0
  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Testing function ABRAM0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram0_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = abram0 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test02

!*****************************************************************************80
!
!! TEST02 tests ABRAM1.
!
  implicit none
!
  real ( kind = 8 ) abram1
  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Testing function ABRAM1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram1_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = abram1 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test03

!*****************************************************************************80
!
!! TEST03 tests ABRAM2.
!
  implicit none

  real ( kind = 8 ) abram2
  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Testing function ABRAM2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram2_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = abram2 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test04

!*****************************************************************************80
!
!! TEST04 tests AIRY_AI_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) airy_ai_int
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Testing function AIRY_AI_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = airy_ai_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test05

!*****************************************************************************80
!
!! TEST05 tests AIRY_BI_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) airy_bi_int
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Testing function AIRY_BI_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = airy_bi_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test06

!*****************************************************************************80
!
!! TEST06 tests AIRY_GI.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) airy_gi
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Testing function AIRY_GI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_gi_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = airy_gi ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test07

!*****************************************************************************80
!
!! TEST07 tests AIRY_HI.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) airy_hi
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Testing function AIRY_HI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_hi_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = airy_hi ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test08

!*****************************************************************************80
!
!! TEST08 tests ARCTAN_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) arctan_int
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Testing function ARCTAN_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctan_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = arctan_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test09

!*****************************************************************************80
!
!! TEST09 tests BESSEL_I0_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bessel_i0_int
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Testing function BESSEL_I0_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = bessel_i0_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test10

!*****************************************************************************80
!
!! TEST10 tests BESSEL_J0_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bessel_j0_int
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Testing function BESSEL_J0_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = bessel_j0_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test11

!*****************************************************************************80
!
!! TEST11 tests BESSEL_K0_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bessel_k0_int
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Testing function BESSEL_K0_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k0_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = bessel_k0_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test12

!*****************************************************************************80
!
!! TEST12 tests BESSEL_Y0_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) bessel_y0_int
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  Testing function BESSEL_Y0_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = bessel_y0_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test13

!*****************************************************************************80
!
!! TEST13 tests CLAUSEN.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) clausen
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Testing function CLAUSEN'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call clausen_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = clausen ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test14

!*****************************************************************************80
!
!! TEST14 tests DEBYE1.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) debye1
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  Testing function DEBYE1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye1_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = debye1 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test15

!*****************************************************************************80
!
!! TEST15 tests DEBYE2.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) debye2
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  Testing function DEBYE2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye2_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = debye2 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test16

!*****************************************************************************80
!
!! TEST16 tests DEBYE3.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) debye3
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  Testing function DEBYE3'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye3_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = debye3 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test17

!*****************************************************************************80
!
!! TEST17 tests DEBYE4.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) debye4
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  Testing function DEBYE4'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye4_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = debye4 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test18

!*****************************************************************************80
!
!! TEST18 tests EXP3_INT.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) exp3_int
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  Testing function EXP3_INT'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exp3_int_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = exp3_int ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test19

!*****************************************************************************80
!
!! TEST19 tests GOODWIN.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  real ( kind = 8 ) goodwin
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  Testing function GOODWIN'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call goodwin_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = goodwin ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test20

!*****************************************************************************80
!
!! TEST20 tests I0ML0.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  real ( kind = 8 ) i0ml0
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  Testing function I0ML0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i0ml0_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = i0ml0 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test21

!*****************************************************************************80
!
!! TEST21 tests I1ML1.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  real ( kind = 8 ) i1ml1
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  Testing function I1ML1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i1ml1_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = i1ml1 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test22

!*****************************************************************************80
!
!! TEST22 tests LOBACHEVSKY.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  real ( kind = 8 ) lobachevsky
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  Testing function LOBACHEVSKY'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lobachevsky_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = lobachevsky ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test23

!*****************************************************************************80
!
!! TEST23 tests STROMGEN.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) stromgen
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST23'
  write ( *, '(a)' ) '  Testing function STROMGEN'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call stromgen_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = stromgen ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test24

!*****************************************************************************80
!
!! TEST24 tests STRUVE_H0.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) struve_h0
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  Testing function STRUVE_H0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_h0_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = struve_h0 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test25

!*****************************************************************************80
!
!! TEST25 tests STRUVE_H1.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) struve_h1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  Testing function STRUVE_H1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_h1_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = struve_h1 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test26

!*****************************************************************************80
!
!! TEST26 tests STRUVE_L0.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) struve_l0
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  Testing function STRUVE_L0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_l0_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = struve_l0 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test27

!*****************************************************************************80
!
!! TEST27 tests STRUVE_L1.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) struve_l1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  Testing function STRUVE_L1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_l1_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = struve_l1 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test28

!*****************************************************************************80
!
!! TEST28 tests SYNCH1.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) synch1
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST28'
  write ( *, '(a)' ) '  Testing function SYNCH1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call synch1_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = synch1 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test29

!*****************************************************************************80
!
!! TEST29 tests SYNCH2.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) synch2
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST29'
  write ( *, '(a)' ) '  Testing function SYNCH2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call synch2_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = synch2 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test30

!*****************************************************************************80
!
!! TEST30 tests TRAN02.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran02
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST30'
  write ( *, '(a)' ) '  Testing function TRAN02'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran02_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran02 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test31

!*****************************************************************************80
!
!! TEST31 tests TRAN03.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran03
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST31'
  write ( *, '(a)' ) '  Testing function TRAN03'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran03_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran03 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test32

!*****************************************************************************80
!
!! TEST32 tests TRAN04.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran04
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST32'
  write ( *, '(a)' ) '  Testing function TRAN04'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran04_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran04 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test33

!*****************************************************************************80
!
!! TEST33 tests TRAN05.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran05
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST33'
  write ( *, '(a)' ) '  Testing function TRAN05'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran05_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran05 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test34

!*****************************************************************************80
!
!! TEST34 tests TRAN06.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran06
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST34'
  write ( *, '(a)' ) '  Testing function TRAN06'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran06_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran06 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test35

!*****************************************************************************80
!
!! TEST35 tests TRAN07.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran07
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST35'
  write ( *, '(a)' ) '  Testing function TRAN07'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran07_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran07 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test36

!*****************************************************************************80
!
!! TEST36 tests TRAN08.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran08
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST36'
  write ( *, '(a)' ) '  Testing function TRAN08'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran08_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran08 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
subroutine test37

!*****************************************************************************80
!
!! TEST37 tests TRAN09.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) comp
  real ( kind = 8 ) fx
  integer n_data
  real ( kind = 8 ) relerr
  real ( kind = 8 ) tran09
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST37'
  write ( *, '(a)' ) '  Testing function TRAN09'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Argument             Abs. error             Rel. error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran09_values ( n_data, x, fx )

    if ( n_data <= 0 ) then
      exit
    end if

    comp = tran09 ( x )
    abserr = abs ( fx - comp )
    relerr = abserr / abs ( fx )
    write ( *, '(2x,f15.10,2x,d15.5,8x,d15.5)' ) x, abserr, relerr

  end do

  return
end
