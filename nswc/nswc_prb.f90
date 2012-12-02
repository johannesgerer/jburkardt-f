program main

!*****************************************************************************80
!
!! MAIN is the main program for NSWC_PRB.
!
!  Discussion:
!
!    NSWC_PRB tests the NSWC library.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NSWC_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the NSWC library.'
 
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
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NSWC_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests AI and BI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ai
  real ap
  real ap2
  real ax
  real ax2
  real bi
  real bp
  real bp2
  real bx
  real bx2
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  AI evaluates the Airy AI function.'
  write ( *, '(a)' ) '  BI evaluates the Airy BI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       AI(X)'
  write ( *, '(a)') ' '
  n = 0

  do

    call airy_values ( n, x, ax, ap, bx, bp )

    if ( n == 0 ) then
      exit
    end if

    ax2 = ai ( x )

    write ( *, '(f8.4,2g14.6)' ) x, ax, ax2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       BI(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call airy_values ( n, x, ax, ap, bx, bp )

    if ( n == 0 ) then
      exit
    end if

    bx2 = bi ( x )

    write ( *, '(f8.4,2g14.6)' ) x, bx, bx2

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BESI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real, parameter :: alpha = 0.0E+00
  real fx
  real fx2(1)
  integer n
  integer nz
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  BESI evaluates the Bessel I function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	   Exact F	 BESI(0)(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call besi0_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0E+00 ) then
      cycle
    end if

    call besi ( x, alpha, 1, 1, fx2(1), nz )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2(1)

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BESI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real, parameter :: alpha = 1.0E+00
  real fx
  real fx2(1)
  integer n
  integer nz
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  BESI evaluates the Bessel I function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       BESI(1)(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call besi1_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0E+00 ) then
      cycle
    end if

    call besi ( x, alpha, 1, 1, fx2(1), nz )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2(1)

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests BESI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real alpha
  real fx
  real fx2(1)
  integer n
  integer nu
  integer nz
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  BESI evaluates the Bessel I function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NU     X       Exact F       BESI(NU)(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call besin_values ( n, nu, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0E+00 ) then
      cycle
    end if

    alpha = real ( nu )
    call besi ( x, alpha, 1, 1, fx2(1), nz )

    write ( *, '(i4,f8.4,2g14.6)' ) nu, x, fx, fx2(1)

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests BESJ.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real, parameter :: alpha = 0.0E+00
  real fx
  real fx2(1)
  integer n
  integer nz
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  BESJ evaluates the Bessel J function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	   Exact F	 BESJ(0)(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call besj0_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0E+00 ) then
      cycle
    end if

    call besj ( x, alpha, 1, fx2(1), nz )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2(1)

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests BESJ.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real, parameter :: alpha = 1.0E+00
  real fx
  real fx2(1)
  integer n
  integer nz
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  BESJ evaluates the Bessel J function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       BESJ(1)(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call besj1_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0E+00 ) then
      cycle
    end if

    call besj ( x, alpha, 1, fx2(1), nz )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2(1)

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests BESJ.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real alpha
  real fx
  real fx2(1)
  integer n
  integer nu
  integer nz
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  BESJ evaluates the Bessel J function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  NU     X       Exact F       BESJ(NU)(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call besjn_values ( n, nu, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x < 0.0E+00 ) then
      cycle
    end if

    alpha = real ( nu )
    call besj ( x, alpha, 1, fx2(1), nz )

    write ( *, '(i4,f8.4,2g14.6)' ) nu, x, fx, fx2(1)

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests BETA.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real beta
  real fxy
  real fxy2
  integer n
  real x
  real y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08:'
  write ( *, '(a)' ) '  BETA evaluates the Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	  Y	   Exact F	 BETA(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call beta_values ( n, x, y, fxy )

    if ( n == 0 ) then
      exit
    end if

    fxy2 = beta ( x, y )

    write ( *, '(2f8.4,2g14.6)' ) x, y, fxy, fxy2

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests BRATIO.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a
  real b
  real fx
  real fx2
  integer ierror
  integer n
  real p
  real q
  real x
  real y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09:'
  write ( *, '(a)' ) '  BRATIO evaluates the normalized incomplete Beta'
  write ( *, '(a)' ) '  function BETAI(A,B,X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A      B       X       Exact F       BETAI(A,B,X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call beta_inc_values ( n, a, b, x, fx )

    if ( n == 0 ) then
      exit
    end if

    y = 1.0E+00 - x

    call bratio ( a, b, x, y, p, q, ierror )

    fx2 = p

    write ( *, '(3f8.4,2g14.6)' ) a, b, x, fx, fx2

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests CIN.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real cin
  real fx
  real fx2
  integer n
  real si
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10:'
  write ( *, '(a)' ) '  CIN evaluates the cosine integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       CIN(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call cin_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    fx2 = cin ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests DAWSON.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real dawson
  real fx
  real fx2
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11:'
  write ( *, '(a)' ) '  DAWSON evaluates Dawson''s integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	   Exact F	 DAWSON(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call dawson_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    fx2 = dawson ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests DILOGARITHM.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real dilogarithm
  real fx
  real fx2
  integer n
  real x
  real x2
  real x3
  real x4
  real fx3
  real fx4

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12:'
  write ( *, '(a)' ) '  DILOGARITHM evaluates the DILOGARITHM function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	   Exact F	 DILOGARITHM(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call dilogarithm_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    fx2 = dilogarithm ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests ELLPF.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a
  real cn
  real dn
  integer ierror
  real k
  real l
  real fx
  real fx2
  integer n
  real sn
  real u
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13:'
  write ( *, '(a)' ) '  ELLPF evaluates the Jacobi elliptic function SN.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A     X       Exact F	    SN(A,X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call jacobi_sn_values ( n, a, x, fx )

    if ( n == 0 ) then
      exit
    end if

    u = x
    k = a
    l = sqrt ( 1.0E+00 - k**2 )

    call ellpf ( u, k, l, sn, cn, dn, ierror )

    fx2 = sn

    write ( *, '(2f8.4,2g14.6)' ) a, x, fx, fx2

  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests ERF.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real erf
  real fx
  real fx2
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14:'
  write ( *, '(a)' ) '  ERF evaluates the ERF function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       ERF(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call erf_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x <= 0.0E+00 ) then
      cycle
    end if

    fx2 = erf ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests EXPLI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  integer ier
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15:'
  write ( *, '(a)' ) '  EXPLI evaluates the exponential integral function EI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       EI(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call ei_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    call expli ( 1, x, fx2, ier )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests EXPLI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  integer ier
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16:'
  write ( *, '(a)' ) '  EXPLI evaluates the exponential integral function E1(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	   Exact F	 E1(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call e1_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    call expli ( 2, x, fx2, ier )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests FRESNEL.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ci
  real fx
  real fx2
  integer n
  real si
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17:'
  write ( *, '(a)' ) '  FRESNEL evaluates the Fresnel cosine integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       C(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call fresnel_cos_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    call fresnel ( x, ci, si )

    fx2 = ci

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FRESNEL evaluates the Fresnel sine integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X	   Exact F	 S(X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call fresnel_sin_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    call fresnel ( x, ci, si )

    fx2 = si

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do
  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests GAMMA.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  real gamma
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18:'
  write ( *, '(a)' ) '  GAMMA evaluates the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       GAMMA(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call gamma_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x <= 0.0E+00 ) then
      cycle
    end if

    fx2 = gamma ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests GAMMP.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a
  real fx
  real fx2
  integer ind
  integer n
  real q
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19:'
  write ( *, '(a)' ) '  GAMMP evaluates the normalized incomplete Gamma'
  write ( *, '(a)' ) '  function P(A,X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A     X       Exact F       GAMMP(A,X)'
  write ( *, '(a)' ) ' '
  n = 0

  do

    call gamma_inc_values ( n, a, x, fx )

    if ( n == 0 ) then
      exit
    end if

    call gratio ( a, x, fx2, q, ind )

    write ( *, '(2f8.4,2g14.6)' ) a, x, fx, fx2

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests MEXP.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 2
  integer, parameter :: lda = n

  real, dimension ( lda, n ) :: a = reshape ( &
    (/ 2.0E+00, 0.0E+00, 2.0E+00, 2.0E+00 /), (/ 2, 2 /) )
  integer ierror
  real, dimension ( lda, n ) :: z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  MEXP computes the matrix exponential.'

  call rmat_print ( lda, n, n, a, '  The matrix A:' )

  call mexp ( a, lda, n, z, ierror )

  call rmat_print ( lda, n, n, z, '  The matrix exponential exp(A):' )

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests PSI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  real psi
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21:'
  write ( *, '(a)' ) '  PSI evaluates the PSI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       PSI(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call psi_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x <= 0.0E+00 ) then
      cycle
    end if

    fx2 = psi ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests SI.
!
!  Modified:
!
!    25 April 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  integer n
  real si
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22:'
  write ( *, '(a)' ) '  SI evaluates the sine integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       SI(X)'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call si_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    fx2 = si ( x )

    write ( *, '(f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
