program main

!*****************************************************************************80
!
!! MAIN is the main program for AIRY_PRB.
!
!  Discussion:
!
!    AIRY_PRB tests some of SLATEC's Airy functions.
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

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the SLATEC library.'

  call ai_prb ( )
  call bi_prb ( )
  call dai_prb ( )
  call dbi_prb ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine ai_prb ( )

!*****************************************************************************80
!
!! AI_PRB demonstrates the use of AI.
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

  real ( kind = 4 ) ai
  real ( kind = 8 ) fx
  real ( kind = 4 ) fx2
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AI_PRB:'
  write ( *, '(a)' ) '  AI computes the Airy function Ai(X)'
  write ( *, '(a)' ) '  using single precision real arithmetic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  AIRY_AI_VALUES returns tabulated values of '
  write ( *, '(a)' ) '  the Airy function Ai(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Ai                       Ai'
  write ( *, '(a)' ) '                      (tabulated)              (computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x, kind = 4 )
    fx2 = ai ( x2 )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine bi_prb ( )

!*****************************************************************************80
!
!! BI_PRB demonstrates the use of BI.
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

  real ( kind = 4 ) bi
  real ( kind = 8 ) fx
  real ( kind = 4 ) fx2
  integer n_data
  real ( kind = 8 ) x
  real ( kind = 4 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BI_PRB:'
  write ( *, '(a)' ) '  BI computes the Airy function Bi(X)'
  write ( *, '(a)' ) '  using single precision real arithmetic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  AIRY_BI_VALUES returns tabulated values of '
  write ( *, '(a)' ) '  the Airy function Bi(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Bi                       Bi'
  write ( *, '(a)' ) '                      (tabulated)              (computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = real ( x, kind = 4 )
    fx2 = bi ( x2 )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine dai_prb ( )

!*****************************************************************************80
!
!! DAI_PRB demonstrates the use of DAI.
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

  real ( kind = 8 ) dai
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DAI_PRB:'
  write ( *, '(a)' ) '  DAI computes the Airy function Ai(X)'
  write ( *, '(a)' ) '  using double precision real arithmetic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  AIRY_AI_VALUES returns tabulated values of '
  write ( *, '(a)' ) '  the Airy function Ai(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Ai                       Ai'
  write ( *, '(a)' ) '                      (tabulated)              (computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dai ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine dbi_prb ( )

!*****************************************************************************80
!
!! DBI_PRB demonstrates the use of DBI.
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

  real ( kind = 8 ) dbi
  real ( kind = 8 ) fx
  real ( kind = 8 ) fx2
  integer n_data
  real ( kind = 8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DBI_PRB:'
  write ( *, '(a)' ) '  DBI computes the Airy function Bi(X)'
  write ( *, '(a)' ) '  using double precision real arithmetic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  AIRY_BI_VALUES returns tabulated values of '
  write ( *, '(a)' ) '  the Airy function Bi(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Bi                       Bi'
  write ( *, '(a)' ) '                      (tabulated)              (computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = dbi ( x )

    write ( *, '(2x,f14.6,2x,g24.16,2x,g24.16)' ) x, fx, fx2

  end do

  return
end
subroutine airy_ai_values ( n_data, x, ai )

!*****************************************************************************80
!
!! AIRY_AI_VALUES returns some values of the Airy Ai(x) function.
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryAi[x]
!
!  Modified:
!
!    11 August 2004
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
!    Output, real ( kind = 8 ) AI, the value of the Airy AI function.
!
  implicit none

  integer, parameter :: n_max = 11

  real ( kind = 8 ) ai
  real ( kind = 8 ), save, dimension ( n_max ) :: ai_vec = (/ &
    0.3550280538878172D+00, &
    0.3292031299435381D+00, &
    0.3037031542863820D+00, &
    0.2788064819550049D+00, &
    0.2547423542956763D+00, &
    0.2316936064808335D+00, &
    0.2098000616663795D+00, &
    0.1891624003981501D+00, &
    0.1698463174443649D+00, &
    0.1518868036405444D+00, &
    0.1352924163128814D+00 /)
  integer n_data
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
    1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    ai = 0.0D+00
  else
    x = x_vec(n_data)
    ai = ai_vec(n_data)
  end if

  return
end
subroutine airy_bi_values ( n_data, x, bi )

!*****************************************************************************80
!
!! AIRY_BI_VALUES returns some values of the Airy Bi(x) function.
!
!  Discussion:
!
!    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
!    solutions of the differential equation:
!
!      W'' - X * W = 0
!
!    In Mathematica, the function can be evaluated by:
!
!      AiryBi[x]
!
!  Modified:
!
!    11 August 2004
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
!    Output, real ( kind = 8 ) BI, the value of the Airy BI function.
!
  implicit none

  integer, parameter :: n_max = 11

  real ( kind = 8 ) bi
  real ( kind = 8 ), save, dimension ( n_max ) :: bi_vec = (/ &
    0.6149266274460007D+00, & 
    0.6598616901941892D+00, & 
    0.7054642029186612D+00, & 
    0.7524855850873156D+00, & 
    0.8017730000135972D+00, & 
    0.8542770431031555D+00, & 
    0.9110633416949405D+00, & 
    0.9733286558781659D+00, & 
    0.1042422171231561D+01, & 
    0.1119872813134447D+01, & 
    0.1207423594952871D+01 /)
  integer n_data
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
    1.0D+00 /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = 0.0D+00
    bi = 0.0D+00
  else
    x = x_vec(n_data)
    bi = bi_vec(n_data)
  end if

  return
end
