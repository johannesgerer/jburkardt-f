program main

!*****************************************************************************80
!
!! MAIN is the main program for DCDFLIB_PRB.
!
!  Discussion:
!
!    DCDFLIB_PRB calls the DCDFLIB tests.
!
!  Modified:
!
!    17 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DCDFLIB_PRB'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DCDFLIB library.'
 
  call test005 ( )
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
  call test24 ( )
  call test25 ( )
  call test26 ( )
  call test27 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DCDFLIB_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test005 ( )

!*****************************************************************************80
!
!! TEST005 tests BETA_INC and BETA_INC_VALUES.
!
!  Modified:
!
!    20 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST005'
  write ( *, '(a)' ) '  BETA_INC computes the incomplete Beta ratio.'
  write ( *, '(a)' ) '  BETA_INC_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    X         Y         A         B         CDF           CDF'
  write ( *, '(a)' ) &
    '                                           (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    y = 1.0D+00 - x

    call beta_inc ( a, b, x, y, cdf_compute, ccdf_compute, ierror )

    write ( *, '(4f10.6,2g14.6)' ) x, y, a, b, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    X         Y         A         B         1-CDF         CCDF'
  write ( *, '(a)' ) &
    '                                           (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    y = 1.0D+00 - x

    call beta_inc ( a, b, x, y, cdf_compute, ccdf_compute, ierror )

    write ( *, '(4f10.6,2g14.6)' ) x, y, a, b, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests CDFBET.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) bound
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CDFBET computes one missing parameter from the'
  write ( *, '(a)' ) '    BETA CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   BETA_CDF ( (P,Q), (X,Y), A, B )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         X         Y         A         B'
  write ( *, '(a)' ) ' '

  do which = 1, 4

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      x = 0.25D+00
      y = 1.0D+00 - x
      a = 2.0D+00
      b = 3.0D+00
    else if ( which == 2 ) then
      p = 0.261719D+00
      q = 1.0D+00 - p
      x = huge ( x )
      y = huge ( y )
      a = 2.0D+00
      b = 3.0D+00
    else if ( which == 3 ) then
      p = 0.261719D+00
      q = 1.0D+00 - p
      x = 0.25D+00
      y = 1.0D+00 - x
      a = huge ( a )
      b = 3.0D+00
    else if ( which == 4 ) then
      p = 0.261719D+00
      q = 1.0D+00 - p
      x = 0.25D+00
      y = 1.0D+00 - x
      a = 2.0D+00
      b = huge ( b )
    end if

    call cdfbet ( which, p, q, x, y, a, b, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFBET returned STATUS = ', status
      cycle
    end if

    write ( *, '(6f10.6)' ) p, q, x, y, a, b

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CDFBIN.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) ompr
  real ( kind = 8 ) p
  real ( kind = 8 ) pr
  real ( kind = 8 ) q
  real ( kind = 8 ) s
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) xn

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CDFBIN computes one missing parameter from the'
  write ( *, '(a)' ) '    Binomial CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   BINOMIAL_CDF ( (P,Q), S, XN, (PR,OMPR) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         S         XN        PR       OMPR'
  write ( *, '(a)' ) ' '

  do which = 1, 4

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      s = 5.0D+00
      xn = 8.0D+00
      pr = 0.875D+00
      ompr = 1.0D+00 - pr
    else if ( which == 2 ) then
      p = 0.067347D+00
      q = 1.0D+00 - p
      s = huge ( s )
      xn = 8.0D+00
      pr = 0.875D+00
      ompr = 1.0D+00 - pr
    else if ( which == 3 ) then
      p = 0.067347D+00
      q = 1.0D+00 - p
      s = 5.0D+00
      xn = huge ( xn )
      pr = 0.875D+00
      ompr = 1.0D+00 - pr
    else if ( which == 4 ) then
      p = 0.067347D+00
      q = 1.0D+00 - p
      s = 5.0D+00
      xn = 8.0D+00
      pr = huge ( pr )
      ompr = huge ( ompr )
    end if

    call cdfbin ( which, p, q, s, xn, pr, ompr, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFBIN returned STATUS = ', status
      cycle
    end if

    write ( *, '(6f10.6)' ) p, q, s, xn, pr, ompr

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests CDFCHI.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) df
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  CDFCHI computes one missing parameter from the'
  write ( *, '(a)' ) '    Chi Square CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   CHI_CDF ( (P,Q), X, DF )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         X         DF'
  write ( *, '(a)' ) ' '

  do which = 1, 3

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      x = 5.0D+00
      df = 8.0D+00
    else if ( which == 2 ) then
      p = 0.242424D+00
      q = 1.0D+00 - p
      x = huge ( x )
      df = 8.0D+00
    else if ( which == 3 ) then
      p = 0.242424D+00
      q = 1.0D+00 - p
      x = 5.0D+00
      df = huge ( df )
    end if

    call cdfchi ( which, p, q, x, df, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFCHI returned STATUS = ', status
      cycle
    end if

    write ( *, '(4f10.6)' ) p, q, x, df

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests CDFCHN.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) df
  real ( kind = 8 ) p
  real ( kind = 8 ) pnonc
  real ( kind = 8 ) q
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  CDFCHN computes one missing parameter from the'
  write ( *, '(a)' ) '    Chi Square CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   CHI_Noncentral_CDF ( (P,Q), X, DF, PNONC )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         X         DF     PNONC'
  write ( *, '(a)' ) ' '

  do which = 1, 4

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      x = 5.0D+00
      df = 8.0D+00
      pnonc = 0.5D+00
    else if ( which == 2 ) then
      p = 0.211040D+00
      q = 1.0D+00 - p
      x = huge ( x )
      df = 8.0D+00
      pnonc = 0.5D+00
    else if ( which == 3 ) then
      p = 0.211040D+00
      q = 1.0D+00 - p
      x = 5.0D+00
      df = huge ( df )
      pnonc = 0.5D+00
    else if ( which == 4 ) then
      p = 0.211040D+00
      q = 1.0D+00 - p
      x = 5.0D+00
      df = 8.0D+00
      pnonc = huge ( pnonc )
    end if

    call cdfchn ( which, p, q, x, df, pnonc, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFCHN returned STATUS = ', status
      cycle
    end if

    write ( *, '(5f10.6)' ) p, q, x, df, pnonc

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests CDFF.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) dfd
  real ( kind = 8 ) dfn
  real ( kind = 8 ) f
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  CDFF computes one missing parameter from the'
  write ( *, '(a)' ) '    F CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   F_CDF ( (P,Q), F, DFN, DFD )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         F         DFN      DFD'
  write ( *, '(a)' ) ' '

  do which = 1, 4

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      f = 5.0D+00
      dfn = 8.0D+00
      dfd = 3.0D+00
    else if ( which == 2 ) then
      p = 0.893510D+00
      q = 1.0D+00 - p
      f = huge ( f )
      dfn = 8.0D+00
      dfd = 3.0D+00
    else if ( which == 3 ) then
      p = 0.893510D+00
      q = 1.0D+00 - p
      f = 5.0D+00
      dfn = huge ( dfn )
      dfd = 3.0D+00
    else if ( which == 4 ) then
      p = 0.893510D+00
      q = 1.0D+00 - p
      f = 5.0D+00
      dfn = 8.0D+00
      dfd = huge ( dfn )
    end if

    call cdff ( which, p, q, f, dfn, dfd, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFF returned STATUS = ', status
      cycle
    end if

    write ( *, '(5f10.6)' ) p, q, f, dfn, dfd

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests CDFFNC.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) dfd
  real ( kind = 8 ) dfn
  real ( kind = 8 ) f
  real ( kind = 8 ) p
  real ( kind = 8 ) pnonc
  real ( kind = 8 ) q
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  CDFFNC computes one missing parameter from the'
  write ( *, '(a)' ) '    noncentral F CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   F_noncentral_CDF ( (P,Q), F, DFN, DFD, PNONC )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         F         DFN       DFD     PNONC'
  write ( *, '(a)' ) ' '

  do which = 1, 5

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      f = 5.0D+00
      dfn = 8.0D+00
      dfd = 3.0D+00
      pnonc = 17.648016D+00
    else if ( which == 2 ) then
      p = 0.60D+00
      q = 1.0D+00 - p
      f = huge ( f )
      dfn = 8.0D+00
      dfd = 3.0D+00
      pnonc = 17.648016D+00
    else if ( which == 3 ) then
      p = 0.60D+00
      q = 1.0D+00 - p
      f = 5.0D+00
      dfn = huge ( dfn )
      dfd = 3.0D+00
      pnonc = 17.648016D+00
    else if ( which == 4 ) then
      p = 0.60D+00
      q = 1.0D+00 - p
      f = 5.0D+00
      dfn = 8.0D+00
      dfd = huge ( dfd )
      pnonc = 17.648016D+00
    else if ( which == 5 ) then
      p = 0.60D+00
      q = 1.0D+00 - p
      f = 5.0D+00
      dfn = 8.0D+00
      dfd = 3.0D+00
      pnonc = huge ( pnonc )
    end if

    call cdffnc ( which, p, q, f, dfn, dfd, pnonc, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFFNC returned STATUS = ', status
      cycle
    end if

    write ( *, '(6f10.6)' ) p, q, f, dfn, dfd, pnonc

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests CDFGAM.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) scale
  real ( kind = 8 ) shape
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  CDFGAM computes one missing parameter from the'
  write ( *, '(a)' ) '    Gamma CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Gamma_CDF ( (P,Q), X, SHAPE, SCALE )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         X         SHAPE    SCALE'
  write ( *, '(a)' ) ' '

  do which = 1, 4

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      x = 5.0D+00
      shape = 8.0D+00
      scale = 3.0D+00
    else if ( which == 2 ) then
      p = 0.981998D+00
      q = 1.0D+00 - p
      x = huge ( x )
      shape = 8.0D+00
      scale = 3.0D+00
    else if ( which == 3 ) then
      p = 0.981998D+00
      q = 1.0D+00 - p
      x = 5.0D+00
      shape = huge ( shape )
      scale = 3.0D+00
    else if ( which == 4 ) then
      p = 0.981998D+00
      q = 1.0D+00 - p
      x = 5.0D+00
      shape = 8.0D+00
      scale = huge ( scale )
    end if

    call cdfgam ( which, p, q, x, shape, scale, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFGAM returned STATUS = ', status
      cycle
    end if

    write ( *, '(5f10.6)' ) p, q, x, shape, scale

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests CDFNBN.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) f
  real ( kind = 8 ) ompr
  real ( kind = 8 ) p
  real ( kind = 8 ) pr
  real ( kind = 8 ) q
  real ( kind = 8 ) s
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  CDFNBN computes one missing parameter from the'
  write ( *, '(a)' ) '    Negative_Binomial CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Negative_BINOMIAL_CDF ( (P,Q), F, S, (PR,OMPR) )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         F         S         PR       OMPR'
  write ( *, '(a)' ) ' '

  do which = 1, 4

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      f = 3.0D+00
      s = 5.0D+00
      pr = 0.875D+00
      ompr = 1.0D+00 - pr
    else if ( which == 2 ) then
      p = 0.988752D+00
      q = 1.0D+00 - p
      f = huge ( f )
      s = 5.0D+00
      pr = 0.875D+00
      ompr = 1.0D+00 - pr
    else if ( which == 3 ) then
      p = 0.988752D+00
      q = 1.0D+00 - p
      f = 3.0D+00
      s = huge ( s )
      pr = 0.875D+00
      ompr = 1.0D+00 - pr
    else if ( which == 4 ) then
      p = 0.988752D+00
      q = 1.0D+00 - p
      f = 3.0D+00
      s = 5.0D+00
      pr = huge ( pr )
      ompr = huge ( ompr )
    end if

    call cdfnbn ( which, p, q, f, s, pr, ompr, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFNBN returned STATUS = ', status
      cycle
    end if

    write ( *, '(6f10.6)' ) p, q, f, s, pr, ompr

  end do

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests CDFNOR.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) mean
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) sd
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  CDFNOR computes one missing parameter from the'
  write ( *, '(a)' ) '    Normal CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   Normal_CDF ( (P,Q), X, MEAN, SD )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         X         MEAN      SD'
  write ( *, '(a)' ) ' '

  do which = 1, 4

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      x = 3.0D+00
      mean = 5.0D+00
      sd = 0.875D+00
    else if ( which == 2 ) then
      p = 0.011135D+00
      q = 1.0D+00 - p
      x = huge ( x )
      mean = 5.0D+00
      sd = 0.875D+00
    else if ( which == 3 ) then
      p = 0.011135D+00
      q = 1.0D+00 - p
      x = 3.0D+00
      mean = huge ( mean )
      sd = 0.875D+00
    else if ( which == 4 ) then
      p = 0.011135D+00
      q = 1.0D+00 - p
      x = 3.0D+00
      mean = 5.0D+00
      sd = huge ( sd )
    end if

    call cdfnor ( which, p, q, x, mean, sd, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFNOR returned STATUS = ', status
      cycle
    end if

    write ( *, '(5f10.6)' ) p, q, x, mean, sd

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests CDFPOI.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) s
  integer ( kind = 4 ) status
  integer ( kind = 4 ) which
  real ( kind = 8 ) xlam

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  CDFPOI computes one missing parameter from the'
  write ( *, '(a)' ) '    Poisson CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   POISSON_CDF ( (P,Q), S, XLAM )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         S         XLAM'
  write ( *, '(a)' ) ' '

  do which = 1, 3

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      s = 3.0D+00
      xlam = 5.0D+00
    else if ( which == 2 ) then
      p = 0.265026D+00
      q = 1.0D+00 - p
      s = huge ( s )
      xlam = 5.0D+00
    else if ( which == 3 ) then
      p = 0.265026D+00
      q = 1.0D+00 - p
      s = 3.0D+00
      xlam = huge ( xlam )
    end if

    call cdfpoi ( which, p, q, s, xlam, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFPOI returned STATUS = ', status
      cycle
    end if

    write ( *, '(4f10.6)' ) p, q, s, xlam

  end do

  return
end
subroutine test11 ( )

!*****************************************************************************80
!
!! TEST11 tests CDFT.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) bound
  real ( kind = 8 ) df
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  integer ( kind = 4 ) status
  real ( kind = 8 ) t
  integer ( kind = 4 ) which

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  CDFT computes one missing parameter from the'
  write ( *, '(a)' ) '    T CDF:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   T_CDF ( (P,Q), T, DF )'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    P         Q         T         DF'
  write ( *, '(a)' ) ' '

  do which = 1, 3

    if ( which == 1 ) then
      p = huge ( p )
      q = huge ( q )
      t = 3.0D+00
      df = 5.0D+00
    else if ( which == 2 ) then
      p = 0.984950D+00
      q = 1.0D+00 - p
      t = huge ( t )
      df = 5.0D+00
    else if ( which == 3 ) then
      p = 0.984950D+00
      q = 1.0D+00 - p
      t = 3.0D+00
      df = huge ( df )
    end if

    call cdft ( which, p, q, t, df, status, bound )

    if ( status /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  CDFT returned STATUS = ', status
      cycle
    end if

    write ( *, '(4f10.6)' ) p, q, t, df

  end do

  return
end
subroutine test12 ( )

!*****************************************************************************80
!
!! TEST12 tests CUMBET and BETA_INC_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x
  real ( kind = 8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  CUMBET computes the Beta CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  BETA_INC_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    X         Y         A         B         CDF           CDF'
  write ( *, '(a)' ) &
    '                                           (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    y = 1.0D+00 - x

    call cumbet ( x, y, a, b, cdf_compute, ccdf_compute )

    write ( *, '(4f10.6,2g14.6)' ) x, y, a, b, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '    X         Y         A         B         1-CDF         CCDF'
  write ( *, '(a)' ) &
    '                                           (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    y = 1.0D+00 - x

    call cumbet ( x, y, a, b, cdf_compute, ccdf_compute )

    write ( *, '(4f10.6,2g14.6)' ) x, y, a, b, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test13 ( )

!*****************************************************************************80
!
!! TEST13 tests CUMBIN and BINOMIAL_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) ompr
  integer ( kind = 4 ) s
  real ( kind = 8 ) s_double
  real ( kind = 8 ) pr
  integer ( kind = 4 ) x
  real ( kind = 8 ) x_double

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  CUMBIN computes the Binomial CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  BINOMIAL_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   X   S    Pr       CDF           CDF'
  write ( *, '(a)' ) '                    (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_cdf_values ( n_data, x, pr, s, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ompr = 1.0D+00 - pr

    s_double = real ( s, kind = 8 )
    x_double = real ( x, kind = 8 )

    call cumbin ( s_double, x_double, pr, ompr, cdf_compute, ccdf_compute )

    write ( *, '(i4,i4,f10.6,2g14.6)' ) s, x, pr, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   X   S    Pr       1-CDF         CCDF'
  write ( *, '(a)' ) '                    (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_cdf_values ( n_data, x, pr, s, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    ompr = 1.0D+00 - pr

    s_double = real ( s, kind = 8 )
    x_double = real ( x, kind = 8 )

    call cumbin ( s_double, x_double, pr, ompr, cdf_compute, ccdf_compute )

    write ( *, '(i4,i4,f10.6,2g14.6)' ) s, x, pr, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test14 ( )

!*****************************************************************************80
!
!! TEST14 tests CUMCHI and CHI_SQUARE_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) df
  real ( kind = 8 ) df_double
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  CUMCHI computes the chi square CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  CHI_SQUARE_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X       DF    CDF           CDF'
  write ( *, '(a)' ) '                 (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_cdf_values ( n_data, df, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    df_double = real ( df, kind = 8 )

    call cumchi ( x, df_double, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,i4,2g14.6)' ) x, df, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X       DF    1-CDF         CCDF'
  write ( *, '(a)' ) '                 (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_cdf_values ( n_data, df, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    df_double = real ( df, kind = 8 )

    call cumchi ( x, df_double, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,i4,2g14.6)' ) x, df, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test15 ( )

!*****************************************************************************80
!
!! TEST15 tests CUMCHN and CHI_NONCENTRAL_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) df
  real ( kind = 8 ) df_double
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  CUMCHN computes the cumulative density'
  write ( *, '(a)' ) '    function for the noncentral chi-squared '
  write ( *, '(a)' ) '    distribution.'
  write ( *, '(a)' ) '  CHI_NONCENTRAL_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DF    Lambda    X         CDF           CDF'
  write ( *, '(a)' ) '                             (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_noncentral_cdf_values ( n_data, x, lambda, df, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    df_double = real ( df, kind = 8 )

    call cumchn ( x, df_double, lambda, cdf_compute, ccdf_compute )

    write ( *, '(2x,i4,2f10.6,2g14.6)' ) &
      df, lambda, x, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    DF    Lambda    X         1-CDF         CCDF'
  write ( *, '(a)' ) '                             (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_noncentral_cdf_values ( n_data, x, lambda, df, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    df_double = real ( df, kind = 8 )


    call cumchn ( x, df_double, lambda, cdf_compute, ccdf_compute )

    write ( *, '(2x,i4,2f10.6,2g14.6)' ) &
      df, lambda, x, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test16 ( )

!*****************************************************************************80
!
!! TEST16 tests CUMF and F_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) dfd
  real ( kind = 8 ) dfd_double
  integer ( kind = 4 ) dfn
  real ( kind = 8 ) dfn_double
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST16'
  write ( *, '(a)' ) '  CUMF computes the F CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  F_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X      DFN DFD    CDF           CDF'
  write ( *, '(a)' ) '                     (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_cdf_values ( n_data, dfn, dfd, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    dfn_double = real ( dfn, kind = 8 )
    dfd_double = real ( dfd, kind = 8 )

    call cumf ( x, dfn_double, dfd_double, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,2i4,2g14.6)' ) x, dfn, dfd, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X      DFN DFD    1-CDF         CCDF'
  write ( *, '(a)' ) '                     (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_cdf_values ( n_data, dfn, dfd, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    dfn_double = real ( dfn, kind = 8 )
    dfd_double = real ( dfd, kind = 8 )

    call cumf ( x, dfn_double, dfd_double, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,2i4,2g14.6)' ) x, dfn, dfd, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test17 ( )

!*****************************************************************************80
!
!! TEST17 tests CUMFNC and F_NONCENTRAL_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) dfd
  real ( kind = 8 ) dfd_double
  integer ( kind = 4 ) dfn
  real ( kind = 8 ) dfn_double
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST17'
  write ( *, '(a)' ) '  CUMFNC computes the noncentral F CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  F_NONCENTRAL_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X      DFN DFD    LAMBDA    CDF           CDF'
  write ( *, '(a)' ) '                               (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_noncentral_cdf_values ( n_data, dfn, dfd, lambda, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    dfn_double = real ( dfn, kind = 8 )
    dfd_double = real ( dfd, kind = 8 )

    call cumfnc ( x, dfn_double, dfd_double, lambda, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,2i4,f10.6,2g14.6)' ) &
      x, dfn, dfd, lambda, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X      DFN DFD    LAMBDA    1-CDF         CCDF'
  write ( *, '(a)' ) '                               (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_noncentral_cdf_values ( n_data, dfn, dfd, lambda, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    dfn_double = real ( dfn, kind = 8 )
    dfd_double = real ( dfd, kind = 8 )

    call cumfnc ( x, dfn_double, dfd_double, lambda, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,2i4,f10.6,2g14.6)' ) &
      x, dfn, dfd, lambda, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test18 ( )

!*****************************************************************************80
!
!! TEST18 tests CUMGAM and GAMMA_INC_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST18'
  write ( *, '(a)' ) '  CUMGAM computes the Gamma CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  GAMMA_INC_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A         X         CDF           CDF'
  write ( *, '(a)' ) '                        (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    call cumgam ( x, a, cdf_compute, ccdf_compute )

    write ( *, '(2f10.6,2g14.6)' ) a, x, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    A         X         CDF           CDF'
  write ( *, '(a)' ) '                        (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    call cumgam ( x, a, cdf_compute, ccdf_compute )

    write ( *, '(2f10.6,2g14.6)' ) a, x, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test19 ( )

!*****************************************************************************80
!
!! TEST19 tests CUMNBN and NEGATIVE_BINOMIAL_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) f
  real ( kind = 8 ) f_double
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) ompr
  integer ( kind = 4 ) s
  real ( kind = 8 ) s_double
  real ( kind = 8 ) pr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST19'
  write ( *, '(a)' ) '  CUMNBN computes the Negative Binomial CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   F   S    Pr       CDF           CDF'
  write ( *, '(a)' ) '                     (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call negative_binomial_cdf_values ( n_data, f, s, pr, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ompr = 1.0D+00 - pr

    f_double = real ( f, kind = 8 )
    s_double = real ( s, kind = 8 )

    call cumnbn ( f_double, s_double, pr, ompr, cdf_compute, ccdf_compute )

    write ( *, '(i4,i4,f10.6,2g14.6)' ) f, s, pr, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   F   S    Pr       1-CDF         CCDF'
  write ( *, '(a)' ) '                     (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call negative_binomial_cdf_values ( n_data, f, s, pr, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    ompr = 1.0D+00 - pr

    f_double = real ( f, kind = 8 )
    s_double = real ( s, kind = 8 )

    call cumnbn ( f_double, s_double, pr, ompr, cdf_compute, ccdf_compute )

    write ( *, '(i4,i4,f10.6,2g14.6)' ) f, s, pr, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test20 ( )

!*****************************************************************************80
!
!! TEST20 tests CUMNOR and NORMAL_01_CDF_VALUES.
!
!  Modified:
!
!    20 September 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  real ( kind = 8 ) mu
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) sigma
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST20'
  write ( *, '(a)' ) '  CUMNOR computes the Normal CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  NORMAL_01_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X         CDF           CDF'
  write ( *, '(a)' ) '              (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    call cumnor ( x, cdf_compute, ccdf_compute )

    write ( *, '(2x,f10.6,2g24.16)' ) x, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X         1-CDF         CCDF'
  write ( *, '(a)' ) '              (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    call cumnor ( x, cdf_compute, ccdf_compute )

    write ( *, '(2x,f10.6,2g24.16)' ) x, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test21 ( )

!*****************************************************************************80
!
!! TEST21 tests CUMPOI and POISSON_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  real ( kind = 8 ) lambda
  integer ( kind = 4 ) n_data
  integer ( kind = 4 ) x
  real ( kind = 8 ) x_double

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST21'
  write ( *, '(a)' ) '  CUMPOI computes the Poisson CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  POISSON_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X    LAMBDA    CDF           CDF'
  write ( *, '(a)' ) '                   (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call poisson_cdf_values ( n_data, lambda, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    x_double = real ( x, kind = 8 )

    call cumpoi ( x_double, lambda, cdf_compute, ccdf_compute )

    write ( *, '(i6,f10.6,2g14.6)' ) x, lambda, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X    LAMBDA    1-CDF         CCDF'
  write ( *, '(a)' ) '                   (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call poisson_cdf_values ( n_data, lambda, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    x_double = real ( x, kind = 8 )
    ccdf_lookup = 1.0D+00 - cdf_lookup

    call cumpoi ( x_double, lambda, cdf_compute, ccdf_compute )

    write ( *, '(i6,f10.6,2g14.6)' ) x, lambda, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test22 ( )

!*****************************************************************************80
!
!! TEST22 tests CUMT and STUDENT_CDF_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) ccdf_compute
  real ( kind = 8 ) ccdf_lookup
  real ( kind = 8 ) cdf_compute
  real ( kind = 8 ) cdf_lookup
  integer ( kind = 4 ) df
  real ( kind = 8 ) df_double
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST22'
  write ( *, '(a)' ) '  CUMT computes the Student T CDF'
  write ( *, '(a)' ) '    and the complementary CDF.'
  write ( *, '(a)' ) '  STUDENT_CDF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X       DF    CDF           CDF'
  write ( *, '(a)' ) '                 (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_cdf_values ( n_data, df, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    df_double = real ( df, kind = 8 )

    call cumt ( x, df_double, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,i4,2g14.6)' ) x, df, cdf_lookup, cdf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X       DF    1-CDF         CCDF'
  write ( *, '(a)' ) '                 (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_cdf_values ( n_data, df, x, cdf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    ccdf_lookup = 1.0D+00 - cdf_lookup

    df_double = real ( df, kind = 8 )

    call cumt ( x, df_double, cdf_compute, ccdf_compute )

    write ( *, '(f10.6,i4,2g14.6)' ) x, df, ccdf_lookup, ccdf_compute

  end do

  return
end
subroutine test23 ( )

!*****************************************************************************80
!
!! TEST23 tests BETA and GAMMA.
!
!  Modified:
!
!    14 April 2007
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
  write ( *, '(a)' ) 'TEST23'
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
subroutine test24 ( )

!*****************************************************************************80
!
!! TEST24 tests ERROR_F, ERROR_FC and ERF_VALUES.
!
!  Modified:
!
!    17 November 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) error_f
  real ( kind = 8 ) erf_compute
  real ( kind = 8 ) erf_lookup
  real ( kind = 8 ) erfc_compute
  real ( kind = 8 ) erfc_lookup
  real ( kind = 8 ) error_fc
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST24'
  write ( *, '(a)' ) '  ERROR_F computes the error function ERF;'
  write ( *, '(a)' ) '  ERROR_FC the complementary error function ERFC.'
  write ( *, '(a)' ) '  ERF_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X         ERF           ERF'
  write ( *, '(a)' ) '              (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, erf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    erf_compute = error_f ( x )

    write ( *, '(f10.6,2g14.6)' ) x, erf_lookup, erf_compute

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X         ERFC          ERFC'
  write ( *, '(a)' ) '              (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  ind = 0
  n_data = 0

  do

    call erf_values ( n_data, x, erf_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    erfc_lookup = 1.0D+00 - erf_lookup
    erfc_compute = error_fc ( ind, x )

    write ( *, '(f10.6,2g14.6)' ) x, erfc_lookup, erfc_compute

  end do

  return
end
subroutine test25 ( )

!*****************************************************************************80
!
!! TEST25 tests GAMMA and GAMMA_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) gamma
  real ( kind = 8 ) gamma_compute
  real ( kind = 8 ) gamma_lookup
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST25'
  write ( *, '(a)' ) '  GAMMA computes the Gamma function;'
  write ( *, '(a)' ) '  GAMMA_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X         GAMMA         GAMMA'
  write ( *, '(a)' ) '              (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, gamma_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    gamma_compute = gamma ( x )

    write ( *, '(f10.6,2g14.6)' ) x, gamma_lookup, gamma_compute

  end do

  return
end
subroutine test26 ( )

!*****************************************************************************80
!
!! TEST26 tests GAMMA_INC and GAMMA_INC_INV.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ), parameter :: test_num = 10

  real ( kind = 8 ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) x
  real ( kind = 8 ) x0
  real ( kind = 8 ) x2

  a = 3.0D+00
  ind = 1
  x0 = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST26'
  write ( *, '(a)' ) '  GAMMA_INC evaluates the incomplete Gamma ratio;'
  write ( *, '(a)' ) '  GAMMA_INC_INV inverts it.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '    A = ', a
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X             P             Q             Inverse'
  write ( *, '(a)' ) ' '

  do i = 0, test_num

    x = dble ( i ) / dble ( test_num )

    call gamma_inc ( a, x, p, q, ind )

    call gamma_inc_inv ( a, x2, x0, p, q, ierror )

    write ( *,     '(4g14.6)' ) x, p, q, x2

  end do

  return
end
subroutine test27 ( )

!*****************************************************************************80
!
!! TEST27 tests PSI and PSI_VALUES.
!
!  Modified:
!
!    14 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) psi
  real ( kind = 8 ) psi_compute
  real ( kind = 8 ) psi_lookup
  integer ( kind = 4 ) n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST27'
  write ( *, '(a)' ) '  PSI computes the Psi function;'
  write ( *, '(a)' ) '  PSI_VALUES looks up some values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    X         PSI           PSI'
  write ( *, '(a)' ) '              (Lookup)      (Computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, psi_lookup )

    if ( n_data == 0 ) then
      exit
    end if

    psi_compute = psi ( x )

    write ( *, '(f10.6,2g14.6)' ) x, psi_lookup, psi_compute

  end do

  return
end
