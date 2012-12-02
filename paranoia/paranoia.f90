program main

!*****************************************************************************80
!
!! MAIN is the main program for PARANOIA.
!
!  Discussion:
!
!    PARANOIA investigates floating point arithmetic.
!
!  Author:
!
!    William Kahan
!
!  Local Parameters:
!
!    Local, integer FROM, the number of the milestone to return to on restart.
!
!    Local, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    Local, integer MILES, the number of milestone reached so far in testing.
!
!    Local, integer START, a flag to tell whether we are restarting or
!    starting from scratch.
!
  implicit none

  real, parameter :: fp0 = 0.0E+00
  real, parameter :: fp1 = 1.0E+00
  real, parameter :: fp2 = 2.0E+00
  real, parameter :: fp27 = 27.0E+00
  real, parameter :: fp3 = 3.0E+00
  real, parameter :: fp4 = 4.0E+00
  real, parameter :: fp8 = 8.0E+00
  real, parameter :: fp9 = 9.0E+00
  integer from
  real, parameter :: half = 0.5E+00
  integer ieee
  integer miles
  real, parameter :: minus1 = -1.0E+00
  real r1
  real r2
  real r3
  real r4
  integer start
  real sticky
  integer temp

  common /roundn/  r1, r2, r3, r4, sticky

  common /global/ fails, sdefct, defect, flaws, radix, ulppls, &
    ulpmin, precis, w, mulgrd, divgrd, subgrd, a1, onemin

  integer fails, sdefct, defect, flaws
  real radix, ulppls
  real ulpmin, precis, w, mulgrd, divgrd, subgrd, a1, onemin
  logical flagf
!
!        fails  ... number of failures
!        sdefct ... number of serious defects
!        defect ... number of defects
!        flaws  ... number of flaws
!
!    Input, integer RADIX, the computed radix of the machine.
!        ulppls ... one unit in the last place (ulp) of 1+epsilon
!        ulpmin ... one unit in the last place (ulp) of 1-epsilon
!        precis ... computed precision of the machine
!        w      ... radix ** precis
!        mulgrd ... used to test for use of guard digits in multiply
!        divgrd ... used to test for use of guard digits in divide
!        subgrd ... used to test for use of guard digits in add/subtract
!        a1     ...
!        onemin ... one minus ulpmin = ~.99999999999999
!
        common /ovunfl/  c1, h1, minpos, nulps, uflthr, phony0
        real c1, h1, minpos, nulps, uflthr, phony0
!        c1     ... 1/c ~= radix^large_integer
!        h1     ... max (2,radix)
!        minpos ... minimum positive number found by mult./division
!        nulps  ... (small integer) * 1 unit in last place of 1+ ... e9
!        uflthr ... the underflow threshold u0
!        phony0 ... candidate and final phony zero if it exists. z0
!
  from = 0

  call timestamp ( )
!
!  Check for restart
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Is this...'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  0  the first time the program is being run'
  write ( *, '(a)' ) '     or'
  write ( *, '(a)' ) '  1  a restart run of the program after failure?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Please enter 0 or 1:'

  read ( *, * ) start

  flagf = ( start /= 0 )

  inquire ( file = 'tst', exist = flagf )

  if ( flagf ) then

    open ( unit = 3, file = 'tst', form = 'unformatted', status = 'old' )
    rewind ( unit = 3 )

    open ( unit = 4, file = 'log', form = 'unformatted', status = 'old' )
    rewind ( unit = 4 )

  else

    open ( unit = 3, file = 'tst', form = 'unformatted', status = 'new' )
    open ( unit = 4, file = 'log', form = 'unformatted', status = 'new' )

  end if

  if ( start /= 0 ) then

    call check_read ( a1, c1, defect, divgrd, fails, flaws, h1, ieee, &
      miles, minpos, mulgrd, nulps, onemin, phony0, precis, r1, r2, &
      r3, r4, radix, sdefct, sticky, subgrd, uflthr, ulppls, ulpmin, w )

    from = miles

    write ( *, '(a)' ) ' '
    write ( *, '(a,i5,a)' ) '  Restarting from milestone ', from, '.'

    if ( from ==   7) go to   881
    if ( from ==  79) go to  3959
    if ( from ==  90) go to  3960
    if ( 105 <= from .and. from <= 109) go to 10100
    if ( from == 115) go to 10100
    if ( 120 <= from .and. from <= 125) go to 10100
    if ( from == 131) go to 10100
    if ( from == 161) go to 10200
    if ( 210 <= from .and. from <= 205) go to 10200
    if ( from == 211) go to 10300
    if ( from == 212) go to 10300
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARANOIA - Error!'
    write ( *, '(a,i6)' ) '  Unrecognized restart milestone "FROM" = ', from
    stop

  end if

10000  continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARANOIA'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  A "paranoid" program to diagnose floating arithmetic.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Using single precision arithmetic.'

  miles = 0
!
!  Count failures, serious defects, defects, flaws
!
  fails = 0
  sdefct = 0
  defect = 0
  flaws = 0
!
!  Print the introduction.
!
  call intro ( )
!
!  Small integer testing.
!
  miles = 7
  call check_write ( ieee, miles )
  miles = miles + 1

881 continue

  call small_int ( defect, fails, fp1, fp2, from, ieee, miles, sdefct )

  from = 0
!
!  Find radix and precision.
!
  call radx ( defect, fails, flaws, fp1, fp2, fp4, fp9, half, &
    miles, onemin, precis, radix, sdefct, ulpmin, ulppls, w )
!
!  Test for extra precision in subexpressions.
!
1680 continue

  miles = 30

  call extra ( fails, fp1, fp2, fp3, fp4, fp8, half, &
    onemin, radix, sdefct, ulpmin, ulppls )

  call check_write ( ieee, miles )
  miles = miles + 1
!
!  Check for guard digits and normalization in subtraction.
!
  miles = 35

  call guard ( defect, divgrd, fails, fp1, fp2, fp3, fp9, &
    half, mulgrd, onemin, radix, sdefct, subgrd, ulpmin, ulppls, w )

  miles = 40
  call check_write ( ieee, miles )
  miles = miles + 1
!
!  Test rounding in multiply, divide, and add/subtract.
!
  call round ( a1, divgrd, fails, flaws, fp1, fp2, fp3, fp9, half, &
    miles, minus1, mulgrd, onemin, r1, r2, r3, radix, sticky, subgrd, &
    ulpmin, ulppls )

  miles = 60
!
!  Test for commutative multiplication.
!
  call commute ( defect, fp3, half, ulpmin, ulppls )
  miles = 70
!
!  Test square root.
!
3959    continue

  call square ( a1, defect, fails, fp1, fp2, fp4, fp8, fp9, from, &
    half, ieee, miles, minus1, onemin, precis, r4, radix, sdefct, &
    ulpmin, ulppls, w )

  from = 0

3960    continue

  miles = 90
  call check_write ( ieee, miles )
  miles = miles + 1
!
!  Test Y to power X.
!
  call power ( a1, defect, from, ieee, miles, w )
  from = 0
!
!  Test underflow thresholds.
!
10100   continue

  call underflow ( c1, fp1, fp2, fp3, fp8, from, h1, half, ieee, &
    miles, minpos, nulps, phony0, uflthr )

  from = 0
  call check_write ( ieee, miles )
  miles = miles + 1
!
!  Test overflow thresholds.
!
10200   continue

  call overflow ( c1, defect, fails, flaws, fp1, fp2, from, h1, &
    half, ieee, miles, minpos, nulps, onemin, phony0, precis, &
    radix, sdefct, uflthr, ulpmin, ulppls, w )

  from = 0
  miles = 210
!
!  See what happens when you divide by zero.
!
10300   continue

  call zeros ( fp0, from, ieee, miles )
  from = 0
  call check_write ( ieee, miles )
  miles = miles + 1

  write ( *, '(a,i4)' ) '  The number of  failures =       ', fails
  write ( *, '(a,i4)' ) '  The number of serious defects = ', sdefct
  write ( *, '(a,i4)' ) '  The number of defects =         ', defect
  write ( *, '(a,i4)' ) '  The number of flaws =           ', flaws

  write ( *, '(a)' ) ' '

  if ( 0 < fails + sdefct + defect + flaws ) then
    go to 6270
  end if

  if ( r1 + r2 + r3 + r4 < fp4 ) then
    go to 6260
  end if

  if ( sticky < fp1 .or. (radix-fp2)*(radix-fp9-fp1) /= fp0 ) then
    go to 6250
  end if

  if ( radix == fp2 .and. &
    ( precis - fp4*fp3*fp2) * (precis - fp27-fp27+fp1) == fp0 ) then
     temp = 754
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Rounding appears to conform to the proposed'
    write ( *, '(a)' ) '  IEEE standard of page 754.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Rounding appears to conform to the proposed'
    write ( *, '(a)' ) '  IEEE standard of page 854.'
  end if

  if ( ieee == 0 ) then
    write ( *, '(a)' ) '  except possibly for single rounding during '
    write ( *, '(a)' ) '  gradual underflow.'
  end if

6250  continue  
    write ( *, '(a)' ) ' '
    write ( *, '(a)' )'  The arithmetic appears to be excellent!'
    go to 6310

6260 continue   
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The arithmetic seems satisfactory.'
    go to 6310

6270    continue

  if ( fails + sdefct + defect == 0 .and. 0 < flaws ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The arithmetic seems satisfactory though flawed.'
  end if

6280    continue

  if ( fails + sdefct == 0 .and. 0 < defect ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The arithmetic may be acceptable despite'
    write ( *, '(a)' ) '  inconvenient defects.'
  end if

6290    continue

  if ( 0 < fails + sdefct ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The arithmetic has unacceptable serious defects.'
  end if

6300    continue

  if ( 0 < fails ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Potentially fatal failure may have spoiled this'
    write ( *, '(a)' ) '  program''s subsequent diagnoses.'
  end if

6310    continue

  close ( unit = 3, status = 'delete' )
  close ( unit = 4, status = 'delete' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARANOIA:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine badsqr ( sflag, z, y )

!*****************************************************************************80
!
!! BADSQR reports on errors involving square of the square root.
!
!  Modified:
!
!    15 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer SFLAG, indicates need to use prefix "SERIOUS".
!
!    Input, real Z, incorrectly computed square of square root.
!
!    Input, real Y, number whose square root squared is wrong.
!
  implicit none

  integer sflag
  real y
  real z

  write ( *, '(a)' ) ' '
  if ( sflag == 1 ) then
    write ( *, '(a)' ) 'Serious defect:'
  else
    write ( *, '(a)' ) 'Defect:'
  end if

  write ( *, '(a,e16.8)' ) '  Comparison alleges that what prints as Z = ', z
  write ( *, '(a,e16.8)' ) '  is too far from sqrt(Z)^2 = ', y

  return
end
subroutine check_read ( a1, c1, defect, divgrd, fails, flaws, h1, ieee, &
  miles, minpos, mulgrd, nulps, onemin, phony0, precis, r1, r2, &
  r3, r4, radix, sdefct, sticky, subgrd, uflthr, ulppls, ulpmin, w )

!*****************************************************************************80
!
!! CHECK_READ reads information from a checkpoint file for restarting.
!
!  Modified:
!
!    15 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Output, integer RADIX, the computed radix of the machine.
!
  implicit none

  real a1
  real c1
  integer defect
  real divgrd
  integer fails
  integer flaws
  real h1
  real half
  integer ieee
  integer miles
  real minpos
  real minus1
  real mulgrd
  real nulps
  real onemin
  real phony0
  real precis
  real r1
  real r2
  real r3
  real r4
  real radix
  integer sdefct
  real sticky
  real subgrd
  real uflthr
  real ulppls
  real ulpmin
  real w

  read ( 4 ) fails, sdefct, defect, flaws, radix, ulppls, ulpmin, precis
  read ( 4 ) w, mulgrd, divgrd, subgrd, a1, onemin
  read ( 4 ) c1, h1, minpos, nulps, uflthr, phony0, ieee
  read ( 4 ) r1, r2, r3, r4, sticky
  read ( 4 ) miles

  rewind 4

  return
end
subroutine check_write ( ieee, miles )

!*****************************************************************************80
!
!! CHECK_WRITE writes information to a checkpoint file for restarting.
!
!  Modified:
!
!    25 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    Input, integer MILES, the number of milestone reached so far in testing.
!
  implicit none

  integer ieee
  integer miles

  common /global/ fails, sdefct, defect, flaws, radix, ulppls, &
              ulpmin, precis, w, mulgrd, divgrd, subgrd, a1, onemin
  integer         fails, sdefct, defect, flaws
  real                              radix, ulppls, &
              ulpmin, precis, w, mulgrd, divgrd, subgrd, a1, onemin

  common /ovunfl/  c1, h1, minpos, nulps, uflthr, phony0
  real c1, h1, minpos, nulps, uflthr, phony0

  common /roundn/  r1, r2, r3, r4, sticky
  real r1, r2, r3, r4, sticky

  write ( 4 ) fails, sdefct, defect, flaws, radix, ulppls, ulpmin, precis
  write ( 4 ) w, mulgrd, divgrd, subgrd, a1, onemin
  write ( 4 ) c1, h1, minpos, nulps, uflthr, phony0, ieee
  write ( 4 ) r1, r2, r3, r4, sticky
  write ( 4 ) miles

  rewind 4

  return
end
subroutine cmpxy ( defect, x, y, z, q, n )

!*****************************************************************************80
!
!! CMPXY compares X and Y=Z**Q for equality.
!
!  Modified:
!
!    15 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
  implicit none

  integer defect
  real, parameter :: fp0 = 0.0E+00
  integer n
  integer q
  real x
  real xx
  real y
  real z

  y = z**q

  if ( y == x ) then

    return

  else if ( z <= fp0 .and. q <= fp0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Warning:'

  else if ( n <= 0 ) then

    defect = defect + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'

  end if

  write ( *, '(a,e16.8,a,i3,a,e16.8)' ) &
    '  Computed (',z, ')^(', q, ') = ', y
  write ( *, '(a,e16.8)' ) '  whereas correct answer is ', x
  xx = y - x
  write ( *, '(a,e16.8)' ) '  Difference is ', xx

  n = n + 1

  return
end
subroutine commute ( defect, fp3, half, ulpmin, ulppls )

!*****************************************************************************80
!
!! COMMUTE tests for commutative multiplication.
!
!  Modified:
!
!    20 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input/output, integer DEFECT, the number of defects.
!
!    Input, real FP3, the value 3.
!
!    Input, real HALF, the value 0.5.
!
!    Input, real ULPMIN, the closest relative separation of real numbers.
!
  implicit none

  integer defect
  integer error
  real, parameter :: fp0 = 0.0E+00
  real fp3
  real half
  integer i
  integer nn
  integer, parameter :: numtry = 20
  real r9
  real ulpmin
  real ulppls
  real x
  real x9
  real y
  real y9
  real z
  real z9

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMMUTE:'
  write ( *, '(a)' ) '  Does multiplication commute?'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4,a)' ) '  Test X*Y=Y*X for ', numtry, ' random pairs.'

  error = 0
  r9 =  sqrt ( fp3 )
  x9 = fp0 / fp3

  do i = 1, numtry

    call random ( x, y, x9, r9 )
    y9 = x9
    call random ( x, y, x9, r9 )
    z = x9 * y9
    y = y9 * x9
    z9 = z - y

    if ( z9 /= fp0 ) then
      error = error + 1
      defect = defect + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Defect:'
      write ( *, '(a,e15.7,a,e15.7)' ) '  X*Y /= Y*X, for X = ', x9, ' Y = ', y9
      write ( *, '(a,e15.7,a,e15.7,a,e15.7)' ) '  X*Y=', z, 'Y*X=', y, &
        'X*Y-Y*X=', z9
      write ( *, '(a,i4)' ) '  Test pair number ', i
    end if

  end do

  x9 = fp0 + half / fp3
  y9 = ( ulppls + ulpmin ) + fp0
  z = x9 * y9
  y = y9 * x9
  z9 = ( fp0 + half / fp3 ) * ( ( ulppls + ulpmin ) + fp0 ) &
    - ( ( ulppls + ulpmin ) + fp0 ) * ( fp0 + half / fp3 )

  if ( z9 /= fp0 ) then
    error = error + 1
    defect = defect + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'
    write ( *, '(a,e15.7,a,e15.7)' ) '  X*Y /= Y*X, for X = ', x9, ' Y = ', y9
    write ( *, '(a,e15.7,a,e15.7,a,e15.7)' ) '  X*Y=', z, 'Y*X=', y, &
      'X*Y-Y*X=', z9
    write ( *, '(a,i4)' ) '  Test pair number ', numtry + 1
  end if

  if ( error == 0 )then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  No commutative failures found in multiplication.'
    return
  end if

  return
end
subroutine extra ( fails, fp1, fp2, fp3, fp4, fp8, half, &
  onemin, radix, sdefct, ulpmin, ulppls )

!*****************************************************************************80
!
!! EXTRA tests for extra precision in subexpressions.
!
!  Modified:
!
!    16 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  integer fails
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  real fp3
  real, parameter :: fp32 = 32.0E+00
  real fp4
  real fp8
  real half
  real onemin
  real q
  real radix
  integer sdefct
  real ulpmin
  real ulppls
  real x
  real x1
  real xx
  real y
  real y1
  real z
  real z1
  real z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXTRA:'
  write ( *, '(a)' ) '  Test for extra precision in subexpressions.'
   
  x =  abs ( ( ( fp4 / fp3 - fp1 ) - fp1 / fp4 ) * fp3 - fp1 / fp4 )

  do

    z2 = x
    x = ( fp1 + ( half * z2 + fp32 * z2 * z2 ) ) - fp1

    if ( z2 <= x .or. x <= fp0 ) then
      exit
    end if

  end do

  y =  abs ( ( fp3 / fp4 - fp2 / fp3 ) * fp3 - fp1 / fp4 )
  z = y
  x = y

  do

    z1 = z
    z = ( fp1 / fp2 - ( ( fp1 / fp2 - ( half * z1 + fp32 * z1 * z1 ) ) &
      + fp1 / fp2 ) ) + fp1 / fp2

    if ( z1 <= z .or. z <= fp0 ) then
      exit
    end if

  end do

  do

    y1 = y
    y = ( half - ( ( half - ( half * y1 + fp32 * y1 * y1 ) ) + half ) ) + half

    if ( y < y1 .and. fp0 < y ) then
      cycle
    end if

    x1 = x
    x = ( ( half * x1 + fp32 * x1 * x1 ) - onemin ) + onemin

    if ( x1 <= x .or. x <= fp0 ) then
      exit
    end if

  end do

  if ( x1 /= y1 .or. x1 /= z1 ) then
    sdefct = sdefct + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, '(a)' ) '  Disagreements among the values X1, Y1, Z1,'
    write ( *, '(a,e15.7,a,e15.7,a,e15.7,a)' ) &
      '  X1 = ', x1, ', Y1 = ', y1, ', Z1 = ', z1, ','
    write ( *, '(a)' ) '  are symptoms of inconsistencies introduced by'
    write ( *, '(a)' ) '  extra-precise evaluation of allegedly "optimized"'
    write ( *, '(a)' ) '  arithmetic subexpressions.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Possibly some part of this test is inconsistent.'
    write ( *, '(a)' ) '  Please notify Karpinski!'

    if ( x1 == ulpmin .or. y1 == ulpmin .or. z1 == ulpmin ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  That feature is not tested further by this program.'
      return
    end if

  end if

  if ( z1 == ulpmin .and. z2 == ulppls ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Subexpressions do not appear to be calculated'
    write ( *, '(a)' ) '  with extra precision.'
    return
  end if

  if ( ulpmin <= z1 .or. ulppls <= z2 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Precision test is inconsistent.'
    write ( *, '(a)' ) '  Please notify Karpinski!'
    xx = z1 - ulpmin
    write ( *, '(a,e15.7)' ) '  ULPMIN =      ', ulpmin
    write ( *, '(a,e15.7)' ) '  Z1 - ULPMIN = ', xx 
    xx = z2 - ulppls
    write ( *, '(a,e15.7)' ) '  ULPPLS =      ', ulppls
    write ( *, '(a,e15.7)' ) '  Z1 - ULPPLS = ', xx 
    return
  end if

  if ( z1 <= fp0 .or. z2 <= fp0 ) then

    write ( *, '(a,f4.0)' )  '  Because of an unusual radix, ', radix
    write ( *, '(a)' ) '  or exact rational arithmetic,'
    write ( *, '(a,e15.7)' ) '  a result Z1 = ', z1
    write ( *, '(a,e15.7)' ) '  or       Z2 = ', z2
    write ( *, '(a)' ) '  of an extra-precision test is inconsistent.'
    write ( *, '(a)' ) '  Please notify Karpinski!'

    if ( z1 == z2 ) then
      write ( *, '(a)' ) '  That feature is not tested further by this program.'
      return
    end if

  end if

  x = z1 / ulpmin
  y = z2 / ulppls
  x = max ( x, y )
  q = -log ( x )
  write ( *, '(a)' )'  Some subexpressions appear to be calculated'
  xx = q / log ( radix )
  write ( *, 1842) xx
  xx = q / log ( fp8 + fp2 )
  write ( *, 1843) xx
1842    format(' extra-precisely, with about   ',e15.7,' extra base b digits,')
1843    format(' that is, roughly ',e15.7,' extra significant decimals.')

  return
end
subroutine guard ( defect, divgrd, fails, fp1, fp2, fp3, fp9, &
  half, mulgrd, onemin, radix, sdefct, subgrd, ulpmin, ulppls, w )

!*****************************************************************************80
!
!! GUARD checks for guard digits and normalization in subtraction.
!
!  Modified:
!
!    18 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  real b9
  integer defect
  real divgrd
  integer fails
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  real, parameter :: fp27 = 27.0E+00
  real fp3
  real fp9
  real half
  real mulgrd
  real onemin
  real r
  real radix
  real s
  integer sdefct
  real subgrd
  real t
  real ulpmin
  real ulppls
  real w
  real x
  real y
  real z

  b9 = radix - ulppls
  mulgrd = 1.0E+00
  divgrd = 1.0E+00
  subgrd = 1.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GUARD'
  write ( *, '(a)' ) '  Check for normalized subtraction,'
  write ( *, '(a)' ) '  and guard digits.'

  if ( fp2 <= radix ) then

    x = w / ( radix * radix )
    y = x + fp1
    z = y - x
    t = z + ulppls
    x = t - z

    if ( x /= ulppls ) then
      fails = fails + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a)' ) '  Subtraction is not normalized properly,'
      write ( *, '(a)' ) '  so X=Y does not imply X+Z=Y+Z!'
    else
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Subtraction appears to be normalized properly.'
    end if

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Checking for guard digits in multiply, divide ' // &
    ' and subtract.'

  y = onemin * fp1
  z = fp1 * onemin
  x = onemin - half
  y = (y-half) - x
  z = (z-half) - x
  x = fp1 + ulppls
  t = x * radix

  r = radix * x
  x = t - radix
  x = x - radix * ulppls
  t = r - radix
  t = t - radix * ulppls
  x = x * ( radix - fp1 )
  t = t * ( radix - fp1 )

  if ( x /= fp0 .or. y /= fp0 .or. z /= fp0 .or. t /= fp0 ) then
    sdefct = sdefct + 1
    mulgrd = fp0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, '(a)' ) '  Multiplication lacks a guard digit,'
    write ( *, '(a)' ) '  violating 1*X = X.'
  end if

  z = radix * ulppls
  x = fp1+z
  y = abs ( ( x + z ) - x * x ) - ulppls
  x = fp1-ulppls
  z = abs ( ( x - ulppls ) - x * x ) - ulpmin

  if ( fp0 < y .or. fp0 < z ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Multiplication gets too many last digits wrong.'
  end if

  y = fp1 - ulppls
  x = fp1 + ulppls
  z = fp1 / y
  y = z - x
  x = fp1 / fp3
  z = fp3 / fp9
  x = x - z
  t = fp9 / fp27
  z = z - t

  if ( x /= fp0 .or. y /= fp0 .or. z /= fp0 ) then
    defect = defect + 1
    divgrd = fp0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'
    write ( *, '(a)' ) '  Division lacks a guard digit.'
    write ( *, '(a)' ) '  so error can exceed 1 ULP,'
    write ( *, '(a)' ) '  or 1/3, 3/9, and 9/27 may disagree.'
  end if

  y = onemin / fp1
  x = onemin - half
  y = ( y - half ) - x
  x = fp1 + ulppls
  t = x / fp1
  x = t - x

  if ( x /= fp0 .or. y /= fp0 ) then
    sdefct = sdefct + 1
    defect = defect - 1 + divgrd
    divgrd = fp0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, '(a)' ) '  Division lacks a guard digit, violating X/1 = X.'
  end if

  x = fp1 / ( fp1 + ulppls )
  y = x - half - half

  if ( fp0 <= y ) then
    sdefct = sdefct + 1
    continue
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Very serious defect:'
    write ( *, '(a)' ) '  Computed value of  1/1.00...001 is not less than 1.'
  end if

  x = fp1 - ulppls
  y = fp1 + radix * ulppls
  z = x * radix
  t = y * radix
  r = z / radix
  s = t / radix
  x = r - x
  y = s - y

  if ( x /= fp0 .or. y /= fp0 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Multiplication and/or division'
    write ( *, '(a)' ) '  gets too many last digits wrong.'
  end if

  y = fp1 - ulpmin
  x = fp1 - onemin
  y = fp1-y
  t = radix - ulppls
  z = radix - b9
  t = radix - t

  if ( x /= ulpmin .or. y /= ulpmin .or. z /= ulppls .or. t /= ulppls ) then
    sdefct = sdefct + 1
    subgrd = fp0
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, '(a)' ) '  Subtraction lacks a guard digit'
    write ( *, '(a)' ) '  so cancellation is obscured.'
  end if

  if ( onemin /= fp1 .and. fp0 <= onemin - fp1 ) then
    sdefct = sdefct + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Very serious defect:'
    write ( *, '(a)' ) '  Comparison alleges  (1-ulpmin) < 1  although'
    write ( *, '(a)' ) '  subtraction yields  (1-ulpmin) - 1 = 0'
    write ( *, '(a)' ) '  thereby vitiating such precautions against'
    write ( *, '(a)' ) '  division by zero as:'
    write ( *, '(a)' ) '   ...  if (x=1.0) then ..... else .../(x-1.0)...'
  end if

  if ( mulgrd * divgrd * subgrd == fp1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  These operations appear to have guard digits'
    write ( *, '(a)' ) '  as they should.'
  end if

  return
end
subroutine intro ( )

!*****************************************************************************80
!
!! INTRO prints the introduction.
!
!  Modified:
!
!    19 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    None
!
  implicit none

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lest this program stop prematurely, i.e. before '
  write ( *, '(a)' ) '  displaying "end of test", try to persuade the '
  write ( *, '(a)' ) '  computer not to terminate execution whenever an error '
  write ( *, '(a)' ) '  such as over/underflow or division by zero occurs, '
  write ( *, '(a)' ) '  but rather to persevere with a surrogate value' 
  write ( *, '(a)' ) '  after, perhaps, displaying some warning. '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  If persuasion avails naught, don''t despair, but run'
  write ( *, '(a)' ) '  this program anyway, to see how many milestones it'
  write ( *, '(a)' ) '  pases.  It should pick up just beyond the last error'
  write ( *, '(a)' ) '  and continue.  If not, it needs further debugging.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Users are invited to help debug and augment this'
  write ( *, '(a)' ) '  program so that it will cope with unanticipated '
  write ( *, '(a)' ) '  and newly found compilers and arithmetic pathologies.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Please send suggestions and interesting results to:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Richard Karpinski'
  write ( *, '(a)' ) '    Computer Center U-76'
  write ( *, '(a)' ) '    University of California'
  write ( *, '(a)' ) '    San Francisco, CA 94143-0704'
  write ( *, '(a)' ) '    USA'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Please include the following information:' 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Precision:   single;'
  write ( *, '(a)' ) '    Version: 31 july 1986;'
  write ( *, '(a)' ) '    Computer:'
  write ( *, '(a)' ) '    Compiler:'
  write ( *, '(a)' ) '    Optimization level:'
  write ( *, '(a)' ) '    Other relevant compiler options:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  BASIC version (c) 1983 by Professor W M Kahan.'
  write ( *, '(a)' ) '  Translated to FORTRAN by T Quarles and G Taylor.'
  write ( *, '(a)' ) '  Modified to ANSI 66/ANSI 77 compatible subset by'
  write ( *, '(a)' ) '  Daniel Feenberg and David Gay.'
  write ( *, '(a)' ) '  You may redistribute this program freely if you'
  write ( *, '(a)' ) '  acknowledge the source.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This program should reveal these characteristics:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  B = radix ( 1, 2, 4, 8, 10, 16, 100, 256, or ... ) .'
  write ( *, '(a)' ) '  P = precision, the number of significant B digits.'
  write ( *, '(a)' ) '  U2 = B/B^P = one ulp (unit in the last place) of'
  write ( *, '(a)' ) '    1.000xxx..'
  write ( *, '(a)' ) '  U1 = 1/B^P = one ulp of numbers a little less than 1.0.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  G1, G2, G3 tell whether adequate guard digits '
  write ( *, '(a)' ) '  are carried, with "1" meaning yes, and "0" no.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    G1 for multiplication,'
  write ( *, '(a)' ) '    G2 for division,'
  write ( *, '(a)' ) '    G3 for subtraction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  R1, R2, R3, R4 tell whether arithmetic is rounded'
  write ( *, '(a)' ) '  or chopped;' 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    0=chopped,'
  write ( *, '(a)' ) '    1=correctly rounded,'
  write ( *, '(a)' ) '   -1=some other rounding;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    R1 for multiplication,'
  write ( *, '(a)' ) '    R2 for division,'
  write ( *, '(a)' ) '    R3 for addition and subtraction,'
  write ( *, '(a)' ) '    R4 for square roots.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S records whether a "sticky bit" is used in rounding.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    0 = no; 1 = yes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  U0 = an underflow threshold.'
  write ( *, '(a)' ) '  E0 and Z0 tell whether underflow is abrupt, gradual '
  write ( *, '(a)' ) '    or fuzzy'
  write ( *, '(a)' ) '  V = an overflow threshold, roughly.'
  write ( *, '(a)' ) '  V0  tells, roughly, whether infinity is represented.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Comparisons are checked for consistency with'
  write ( *, '(a)' ) '  subtraction, and for contamination by pseudo-zeros.'
  write ( *, '(a)' ) '  SQRT is tested.  So is Y^X for (mostly) integers X.'
  write ( *, '(a)' ) '  Extra-precise subexpressions are revealed but not '
  write ( *, '(a)' ) '  yet tested.  Decimal-binary conversion is not yet '
  write ( *, '(a)' ) '  tested for accuracy.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program attempts to discriminate between:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    FLAWS,           like the lack of a sticky bit,'
  write ( *, '(a)' ) '    SERIOUS DEFECTS, like the lack of a guard digit,'
  write ( *, '(a)' ) '    FAILURES,        like 2+2 = 5.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The diagnostic capabilities of this program go beyond'
  write ( *, '(a)' ) '  an earlier program called MACHAR, which can be found'
  write ( *, '(a)' ) '  at the end of the book:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    W J Cody and W Waite,'
  write ( *, '(a)' ) '    Software Manual for the Elementary Functions,'
  write ( *, '(a)' ) '    1980.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Although both programs try to discover the radix (B),'
  write ( *, '(a)' ) '  precision (P) and range (over/underflow thresholds)'
  write ( *, '(a)' ) '  of the arithmetic, this program tries to cope with a'
  write ( *, '(a)' ) '  wider variety of pathologies and to say how well the'
  write ( *, '(a)' ) '  arithmetic is implemented.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The program is based upon a conventional radix'
  write ( *, '(a)' ) '  representation for floating-point numbers,'
  write ( *, '(a)' ) '  but also allows for logarithmic encoding (B = 1)'
  write ( *, '(a)' ) '  as used by certain early Wang machines.'

  return
end
subroutine newd ( x, z1, q, z, d, fp1, fp2, half, radix )

!*****************************************************************************80
!
!! NEWD updates D and Z.
!
!  Discussion:
!
!    This subroutine updates D and Z, setting 
!
!      newd = radix * d
!
!    and
!
!      newz^2 mod newd = z^2 mod d
!
!  Modified:
!
!    14 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  real d
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  real half
  real radix
  real temp
  real temp1
  real q
  real x
  real z
  real z1

  x = z1 * q
  temp = half - x / radix
  temp1 = aint ( temp )

  if ( temp < temp1 ) then
    temp1 = temp1 - fp1
  end if

  x = temp1 * radix + x
  q = ( q - x * z ) / radix + x * x * ( d / radix )
  z = z - fp2 * x * d

  if ( z <= fp0 ) then
    z = -z
    z1 = -z1
  end if

  d = radix * d

  return
end
subroutine overflow ( c1, defect, fails, flaws, fp1, fp2, from, h1, &
  half, ieee, miles, minpos, nulps, onemin, phony0, precis, radix, &
  sdefct, uflthr, ulpmin, ulppls, w )

!*****************************************************************************80
!
!! OVERFLOW tests overflow threshholds.
!
!  Discussion:
!
!    Here are the characteristics of a few "classic" computers:
!
!      CDC        sign-magnitude,  infinity symbol, no trap
!      ELXSI      sign-magnitude,  trap or infinity symbol if no trap
!      IBM        sign-magnitude,  saturation value, no trap
!      PRIME      twos-complement, saturation value, no trap
!      VAX        sign-magnitude,  trap
!
!  Modified:
!
!    15 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input/output, integer FROM, counts the milestones previously reached.
!
!    ?, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    ?, integer MILES, the number of milestones passed.
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  real c1
  integer defect
  integer fails
  integer flaws
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  integer from
  real h1
  real half
  integer i
  integer i0
  integer ieee
  integer miles
  real minpos
  real nulps
  real onemin
  real phony0
  real precis
  real radix
  real sat
  integer sdefct
  real temp
  real uflthr
  real ulpmin
  real ulppls
  real v
  real v9
  real w
  real x
  real y
  real z
  integer zflag
!
!  Reassign values to variables using the log file,
!  then go to restart point
!
  if ( from /= 0 )then
    read(3) i, v, sat, v9, x, y, z, zflag
    rewind 3
    if (from == 161) go to 5582
    if (from == 170) go to 5680
    if (from == 175) go to 5810
    if ( 201 <= from .and. from <= 205) go to 5999
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OVERFLOW - Error!'
    write ( *, '(a)' ) '  Unrecognized restart milestone.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OVERFLOW'
  write ( *, '(a)' ) '  Searching for overflow threshold:'

  miles = 161
  call check_write ( ieee, miles )
!
!  Set Y to -1 * a large power of the radix.
!
  y = -c1
  v9 = h1 * y
!
!  Multiply by radix (h1) until overflow occurs.
!
  do

    v = y
    y = v9
    write(3) i,v,sat,v9,x,y,z,zflag
    rewind 3
    v9 = h1 * y

    if ( y <= v9 ) then
      exit
    end if

  end do
!
!  System does not trap on overflow
!
!  Possibilities:
!
!    y < v9,  
!    v9 is the value returned after overflow
!    y is the largest power of the radix
!    v is the second largest power of the radix
!
!    v9 == y, 
!    both are saturation value
!    v is the largest power of the radix
!
!    v9 == y, 
!    both are infinity symbol
!    v is the largest power of the radix
!
!  Test 1: value returned after overflow shrinks in magnitude
!
5540    continue

  if ( v9 /= y ) then
    sdefct = sdefct + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, 5541) y, v9
5541 format('  overflow past  ', 1pe16.8, &
                '  shrinks to ', 1pe16.8)
  end if
!
!  Test 2: 
!  Two's complement machine saturates at negative largest power of the radix.
!  Need to distinguish system with overflow symbols from one without them.
! 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Can " Z = -Y " overflow?'
  write ( *, '(a,1pe16.8)' ) '  Trying it on Y = ', y

  sat = -y

  if ( v - y /= v + sat ) then
    flaws = flaws + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Flaw:'
    write ( *, '(a)' ) '  -(-Y) differs from Y.'
  else
    write ( *, '(a)' )  '  Seems OK.'
  end if

  go to 5590
!
!  Restart point for systems that trap on overflow
!
!  v9 = y =  -(largest power of radix)
!  v      =  -(second largest power of radix)
!
!  Test 2: two's complement machine
!
5582    continue

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Can Z = -Y overflow?'
  write ( *, '(a,1pe16.8)' ) '  Trying it on Y = ', y
!
!  Put something here to handle the trap.
!
  sat = -y

  if (v - y == v + sat ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Seems OK.'
  else
    flaws = flaws + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Flaw:'
    write ( *, '(a)' ) '  -(-Y) differs from Y.'
  end if
!
!  This code works for a sign-magnitude machine, but fails for a 
!  twos-complement one.
!
  v = y * ( h1 * ulppls - h1 )
  v = v + y * ( ( fp1 - h1 ) * ulppls )

  write ( *, '(a,1pe16.8)' ) '  Overflow threshold is  v = ', v
  write ( *, '(a)' ) '  There is no saturation value'
  write ( *, '(a)' ) '  because the system traps on overflow.'
  go to 5640
!
!  non-trapping systems (continued)
!
5590    continue

  y = v * ( h1 * ulppls - h1 )
  z = y + v * ( ( fp1 - h1 ) * ulppls )
  if ( z < sat ) then
    y = z
  end if
  if ( y < sat ) then
    v = y
  end if
!
!  The overflow threshold equals the saturation value
!  if the latter behaves as a number rather than an
!  overflow symbol.  an overflow symbol is not
!  changed when any number is added to or subtracted
!  from it.
!
  if ( sat - v < sat ) then
    v = sat
  end if

  write ( *, '(a,1pe16.8)' ) '  Overflow threshold is  v = ', v
  write ( *, '(a,1pe16.8)' ) '  Overflow saturates at  sat = ', sat

5640    continue

  miles = 163
  write ( *, '(a)' ) '  No overflow should be signaled for  v*1 = '
  temp = v * fp1
  write ( *, '(a,1pe16.8)' ) '                                           ', temp
  write ( *, '(a)' ) '                            nor for  v/1 = '
  temp = v / fp1
  write ( *, '(a,1pe16.8)' ) '                                           ', temp
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Any overflow signal separating this * '
  write ( *, '(a)' ) '  from one above is a defect.'
!
!  Need to add code here to handle overflows just above
!
  miles = 170
!
!  problem: sat not defined if we trapped on overflow above
!
5660    continue

  if ( v <= -v .or. sat <= -sat .or. v <= -uflthr .or. v <= uflthr ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Comparisons are confused by overflow.'
  end if

5680    continue

  miles = 175
  i = 0
  z = uflthr

5700    continue

  i = i + 1

  if ( z /= fp0 ) then

    v9 = sqrt ( z )
    y = v9 * v9

    if ( y / ( fp1 - radix * nulps ) < z .or. &
        ( fp1 + radix * nulps ) * z < y ) then

      if ( v9 <= ulpmin ) then
        zflag = 0
        defect = defect + 1
      else
        zflag = 1
        sdefct = sdefct + 1
      end if

      call badsqr ( zflag, z, y )

    end if

  end if

  if ( i == 1 ) then
    z = minpos
    go to 5700
  else if ( i == 2 ) then
    z = phony0
    go to 5700
  else
    i = 0
    z = v
  end if

5810    continue

  miles=180
!
!  Failure: attempts to evaluate  sqr(overflow threshold v)  in double
!  precision in BASIC on the  ibm pc  display the word  " overflow "
!  and then disable the keyboard!  This is disastrous.
!
  if ( radix == 2. .and. precis == 56.0E+00 .and. &
      phony0 /= 0.0E+00 .and. -fp0 /= fp0 ) then

    fails = fails + 1

  else

    do

      v9 =  sqrt ( z )
      x = ( fp1 - radix * nulps ) * v9
      v9 = v9 * x

      if ( v9 < ( fp1 - fp2 * radix * nulps ) * z .or. z < v9 ) then

        y = v9

        if ( w <= x ) then
          zflag = 0
          defect = defect + 1
        else
          zflag = 1
          sdefct = sdefct + 1
        end if

        i = 1
        call badsqr ( zflag, z, y ) 

      end if

      if ( i == 1 ) then
        exit
      end if

      i = 1
      z = sat
 
    end do

  end if

  miles = 190
  call check_write ( ieee, miles )
  miles = miles + 1

  x = uflthr * v
  y = radix * radix

  if (x * y < fp1 .or. y < x ) then
    if ( x * y < ulpmin .or. y / ulpmin < x ) then
      defect = defect + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Defect:'
      write ( *, '(a)' ) '  Badly unbalanced range.'
      write ( *, 5961) x
5961  format('  uflthr * v  =', e13.5,' is too far from 1.')
    else
      flaws = flaws + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Flaw:'
      write ( *, '(a)' ) '  Unbalanced range.'
      write ( *, 5971) x
5971  format('  uflthr * v  =', e13.5,' is too far from 1.')
    end if
  end if
!
!  Test  x/x   vs.  1
!
  i0 = 1
  go to 6000

5999  continue

  i0 = from - 200

6000    continue

  do i = i0, 5

    if ( i == 2 ) then
      x = fp1 + ulppls
    else if ( i == 3 ) then
      x = v
    else if ( i == 4 ) then
      x = uflthr
    else if ( i == 5 ) then
      x = radix
    else
      x = onemin
    end if

    miles = 200 + i

    if ( from == miles ) then
      sdefct = sdefct + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Serious defect:'
      write ( *, '(a,e13.5)' ) '  x/x traps when x = ', x
      cycle
    end if

    y = x
    call check_write ( ieee, miles )
    v9 = ( y / x - half ) - half

    if ( v9 == fp0 ) then
      cycle
    end if

    if ( v9 == - ulpmin .and. i < 5 ) then
      sdefct = sdefct + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Serious defect:'
      write ( *, '(a,e13.5)' ) '  X/X differs from 1 when X = ', x
      write ( *, '(a,e13.5)' ) '  Instead,  X/X - 1/2 - 1/2 = ', v9

    else
      fails = fails + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a,e13.5)' ) '  X/X differs from 1 when X = ', x
      write ( *, '(a,e13.5)' ) '  Instead,  X/X - 1/2 - 1/2 = ', v9
    end if

  end do

  return
end
subroutine partuf ( z, zname, defect, fp1, fp2, ieee, miles, partu, &
  radix, restrt, sdefct, ulppls )

!*****************************************************************************80
!
!! PARTUF tests for partial underflow.
!
!  Modified:
!
!    23 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    ?, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    ?, integer MILES, the number of milestone reached so far in testing.
!
!    Output, integer PARTU, indicates the presence of partial underflow.
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  integer defect
  real divtmp
  real dummy
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  integer ieee
  integer miles
  real multp1
  real multp2
  integer partu
  real radix
  logical restrt
  integer sdefct
  real ulppls
  real z
  character ( len = 8 ) zname

  partu = 0

  if ( .not. restrt ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PARTUF:'
    write ( *, '(a)' ) '  Testing for partial underflow.'

    call check_write ( ieee, miles )

    if ( z == fp0 ) then
     return
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,a8,a)' ) '  Comparison denies  ', zname, ' = 0.'
    write ( *, '(a)' ) '  Hence, it should be safe to evaluate:'
    write ( *, '(a,a8,a,a8,a,a8)' ) '  ( ', zname, ' + ', zname, ' ) / ', zname
    dummy = ( z + z ) / z
    write ( *, '(a)' ) ' '
    write ( *, '(a,e15.7)' ) '  The computed result is ', dummy

    if ( abs ( dummy - fp2 ) < radix * ulppls ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  This is OK provided over/underflow '
      write ( *, '(a)' ) '  has not just been signaled.'
    else if ( dummy < fp1 .or. fp2 < dummy ) then
      partu = 1
      sdefct = sdefct + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Serious defect!'
    else
      partu = 1
      defect = defect + 1
      write ( *, '(a)' ) '  This is a defect.'
    end if

  else if ( restrt ) then

    partu = 1
    sdefct = sdefct + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Serious defect!'

  end if

  multp1 = z * fp1
  multp2 = fp1 * z
  divtmp = z / fp1

  if ( z /= multp1 .or. z /= multp2 .or. z /= divtmp ) then
    partu = 1
    defect = defect + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'
    write ( *, '(a,a8,a,e16.8,a)' ) &
      '  What prints as  ', zname, ' = ', z, ' compares different from:'
  end if

  if ( z /= multp1 ) then
    write ( *, '(a,a8,a,e16.8)' ) '           ', zname, '*1 = ', multp1
  end if

  if ( z /= multp2 .and. multp2 /= multp1 ) then
    write ( *, '(a,a8,a,e16.8)' ) '           1*', zname, ' = ', multp2
  end if

  if ( z /= divtmp ) then
    write ( *, '(a,a8,a,e16.8)' ) '           ', zname, '/1 = ', divtmp
  end if

  if ( multp2 /= multp1 ) then
    defect = defect + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'
    write ( *, '(a)' ) '  Failure of multiplication commutation.'
    write ( *, '(a)' ) '  The following two results are apparently different:'
    write ( *, '(a,a8,a,e16.8)' ) '  ', zname, '*1 =', multp1
    write ( *, '(a,a8,a,e16.8)' ) '  1*', zname, ' =', multp2
  end if

  if ( 0 < partu ) then
    call check_write ( ieee, miles )
    miles = miles + 1
  end if

  return
end
subroutine power ( a1, defect, from, ieee, miles, w )

!*****************************************************************************80
!
!! POWER tests the calculation of Y**X.
!
!  Modified:
!
!    19 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    ?, integer FROM, the number of the milestone to return to on restart.
!
!    ?, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    ?, integer MILES, the number of milestone reached so far in testing.
!
  implicit none

  real a
  real a1
  integer defect
  real, parameter :: fp0 = 0.0E+00
  real, parameter :: fp1 = 1.0E+00
  real fp2
  real fp3
  real fp4
  real fp8
  integer from
  integer i
  integer ieee
  integer m
  integer miles
  real minus1
  integer n
  integer n1
  integer numtry
  real w
  real x
  real y
  real z

  fp2 = fp1 + fp1
  fp3 = fp2 + fp1
  fp4 = fp3 + fp1
  fp8 = fp4 + fp4
  minus1 = -fp1
  a = fp1 / a1
  numtry = 20

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POWER'
  write ( *, '(a)' ) '  Testing powers Z^I for small integers Z and I.'

  n = 0
  m = 3

  if ( from == 0 ) then

    write ( *, '(a)' ) '  Start with 0.0**0'
    z = -fp0
    i = 0

  else if ( from == 90 ) then

    z = minus1
    i = -4

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'POWER - Error!'
    write ( *, '(a)' ) '  Unrecognized restart milestone.'
    stop

  end if
!
!  Test powers of zero.
!
  do

    x = fp1
    call power_comp ( defect, x, z, i, m, n, w )

    if ( i <= 10 ) then 
      i = 1023
      call power_comp ( defect, x, z, i, m, n, w )
    end if

    if ( z == minus1 ) then
      exit
    end if
!
!  If (-1)^n is invalid, replace 'MINUS1' by 'FP1'
!
    z = minus1
    i = -4

  end do
!
!  Print N if 0 < N.
!
  if ( 0 < n ) then
    write ( *, '(a,i4,a)' ) '  Similar discrepancies have occurred ', n, &
      ' times.'
  end if 

  n1 = n
  n = 0
  z = a1
  m = int ( fp2 * log ( w ) / log ( a1 ) )
!
!  loop
!
  do

    x = z
    i = 1
    call power_comp ( defect, x, z, i, m, n, w )

    if ( z == a ) then
      exit
    end if

    z = a

  end do
!
!  Powers of radix B have been tested; next try a few primes.
!  
  miles = 100
  m = numtry
  z = fp3

  do

    x = z
    i = 1
    call power_comp ( defect, x, z, i, m, n, w )

    do
      z = z + fp2

      if ( fp3 * aint ( z / fp3 ) /= z ) then
        exit
      end if

    end do

    if ( fp8 * fp3 <= z ) then
      exit
    end if

  end do

  if ( 0 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An error like this may invalidate financial'
    write ( *, '(a)' ) '  calculations involving interest rates.'
  end if

  if ( 0 < n ) then
    write ( *, '(a,i4,a)' ) '  Similar discrepancies have occurred ', n, &
      ' times.'
  end if

  n = n + n1

  if ( n == 0 ) then
    write ( *, '(a)' ) '  No discrepancies found.'
  else
    call check_write ( ieee, miles )
    miles = miles + 1
  end if

  return
end
subroutine power_comp ( defect, x, z, i, m, n, w )

!*****************************************************************************80
!
!! POWER_COMP compares Z**I with Z*Z*...*Z.
!
!  Modified:
!
!    19 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
  implicit none

  integer defect
  integer i
  integer m
  integer n
  integer q
  real w
  real x
  real y
  real z

  do

    y = z**i
    q = i
!
!  Test whether Y = X.
!
    call cmpxy ( defect, x, y, z, q, n )

    i = i + 1
!
!  With X = Z^M
!
    if ( m < i ) then
      return
    end if

    x = z * x

    if ( w <= x ) then
      exit
    end if

  end do

  return
end
subroutine radx ( defect, fails, flaws, fp1, fp2, fp4, fp9, half, &
  miles, onemin, precis, radix, sdefct, ulpmin, ulppls, w )

!*****************************************************************************80
!
!! RADX finds the radix and the precision.
!
!  Modified:
!
!    16 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    ?, integer MILES, the number of milestone reached so far in testing.
!
!    Output, integer RADIX, the computed radix of the machine.
!
  implicit none

  real b9
  integer defect
  integer fails
  integer flaws
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  real fp4
  real fp9
  real half
  integer i
  integer miles
  real onemin
  real precis
  real radix
  real radsav
  integer sdefct
  real sixth
  real, parameter :: t8 = 240.0E+00
  real third
  real ulpmin
  real ulpmsv
  real ulppls
  real w
  real x
  real y
  real z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RADX:'
  write ( *, '(a)' ) '  Searching for radix and precision.'
!
!  Looking for W large enough to make 1 insignificant at the
!  precision available.  Increase by powers of 2 until we find it.
!
  w = fp1

  do

    w = w + w
    y = w + fp1
    z = y - w
    y = z - fp1

    if ( ( fp0 <= - fp1 + abs ( y ) ) ) then
      exit
    end if

  end do
!
!  Now W is just big enough that 1 <= |((W+1)-W)-1| ...
!  i.e. 1 is insignificant relative to W.
!
  precis = fp0
  y = fp1

  do

    radix = w + y
    y = y + y
    radix = radix - w

    if ( radix /= fp0 ) then
      exit
    end if

  end do
!
!  If base is 1, it is not characterized by a precision, so
!  don't bother to hunt for it.
!
  if ( radix < fp2 ) then

    radix = fp1
    write ( *, '(a,f7.0)' ) '  Radix = ', radix
!
!  If 2 <= radix, try to find the precision (# significant digits)
!
  else

    write ( *, '(a,f4.0)' ) '  Radix = ', radix
    w = fp1

    do

      precis = precis + fp1
      w = w * radix
      y = w + fp1
      z = y - w

      if ( z /= fp1 ) then
        exit
      end if

    end do

  end if
!
!  Now W = radix^precis is barely too big to satisfy (w+1)-w = 1.
!
  ulpmin = fp1 / w
  ulppls = radix * ulpmin
  write ( *, '(a,1pe16.8)' ) '  Closest relative separation found is ', ulpmin
  write ( *, '(a)' )
  write ( *, '(a)' ) '  Recalculating radix and precision...'
  radsav = radix
  ulpmsv = ulpmin

  x = fp4 / 3.0E+00
  third = x - fp1
  sixth = ( fp1 / fp2 ) - third
  x = sixth + sixth
  x= abs ( x - third )
  x = max ( x, ulppls )
!
!  Now X = (unknown no.) ulps of  1 + ...
!
  do

    ulppls = x
    y = ( fp1 / fp2 ) * ulppls + 32.0E+00 * ulppls * ulppls
    y = fp1 + y
    x = y - fp1
!
!  X => ((x/2) + epsilon) mod (1 ulp of 1+)
!
    if ( ulppls <= x .or. x <= fp0 ) then
      exit
    end if

  end do
!
!  If X does not underflow to 0, then it is still (unknown) * ulp
!  so try again....  otherwise, previous value (ulppls) is 1 ulp
!  ... now  ulppls = 1 ulp of  1 + ...
!
  x = fp2 / 3.0E+00
  sixth = x - ( fp1 / fp2 )
  third = sixth + sixth
  x = third - ( fp1 / fp2 )
  x = abs ( x + sixth )

  x = max ( x, ulpmin )
!
!  Now X = (unknown no.) ulps of  1 - ...
!
  do

    ulpmin = x
    y = ( fp1 / fp2 ) * ulpmin + 32.0E+00 * ulpmin * ulpmin
    y = ( fp1 / fp2 ) - y
    x = ( fp1 / fp2 ) + y
    y = ( fp1 / fp2 ) - x
    x = ( fp1 / fp2 ) + y
!
!  X => (x/2 = epsilon) mod (1 ulp of 1-)
!
    if ( ulpmin <= x .or. x <= fp0 ) then
      exit
    end if

  end do
!
!  Now ULPMIN = 1 ulp of  1 - ...
!
!  Summarize the results.
!
  if ( ulpmin == ulpmsv ) then
    write ( *, '(a)' ) '  ...confirms closest relative separation.'
  else
    write ( *, '(a,e13.5)' ) &
      '  ...gets better closest relative separation = ', ulpmin
  end if

  w = fp1 / ulpmin
  onemin = ( half - ulpmin ) + half
!
!  ... = 1 - ulpmin = nextafter(1.0, 0)
!
  radix = aint ( 0.01E+00 + ulppls / ulpmin )

  if ( radix == radsav ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Radix confirmed.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Mystery:'
    write ( *, '(a,e13.5)' ) '  Recalculated radix = ', radix
  end if
!
!  Radices 1, 2 and 10 pass muster.
!
  if ( radix /= fp2 .and. radix /= 10.0E+00 .and. radix /= fp1 ) then

    if (radix <= 16.0E+00 ) then
      flaws = flaws + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Flaw:'
      write ( *, '(a,f4.0,a)' ) '  Radix ', radix, ' is not as good as 2 or 10.'
    else
      defect = defect + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Defect:'
      write ( *, '(a,f4.0,a)' ) '  Radix ', radix, ' is so big that roundoff'
      write ( *, '(a)' ) 'propagates capriciously.'
    end if

  end if

  miles = 20
!
!  Test for fuzzy comparison.
!
  if ( half <= onemin - half ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  (1-u1)-1/2 < 1/2  is false,'
    write ( *, '(a)' ) '  so this program may malfunction.'
  end if

  x = onemin
  i = 1

  do

    y = x-half
    z = y-half

    if ( x == fp1 .and. z /= fp0 ) then
      fails = fails + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a)' ) '  Comparison is fuzzy; it alleges x=1 although'
      write ( *, '(a,d16.8)' ) '  subtraction yields (x - 1/2) - 1/2 = ', z
    end if

    if ( i == 0 ) then
      exit
    end if

    x = fp1 + ulppls
    i = 0

  end do

  miles = 25
!
!  End of test for fuzzy comparison.
!
  b9 = radix - fp1
  b9 = ( b9 - ulppls ) + fp1

  if ( radix /= fp1 ) then
!
!  b9 = nextafter(radix, 0)
!
    x = - t8 * log ( ulpmin ) / log ( radix )
    y = aint ( half + x )

    if ( abs ( x - y ) * fp4 < fp1 ) then
      x = y
    end if

    precis = x / t8
    y = aint ( half + precis )

    if ( abs ( precis - y ) * t8 < half ) then
      precis = y
    end if
!
!  Purify integers.
!
    if ( precis == aint ( precis ) ) then
      write ( *, '(a,f4.0,a,f6.2)' ) &
        '  The number of significant digits of radix ', radix, ' is ' , precis
      go to 1650
    end if

  end if

  write ( *, '(a)' ) '  Precision cannot be characterized by an integer'
  write ( *, '(a)' ) '  number of significant digits;'
  if ( fp1 < radix ) then
    write ( *, '(a)' ) '  but, by itself, this is a minor flaw.'
    write ( *, '(a,f4.0,a,f6.2)' ) &
      '  The number of significant digits of radix ', radix, ' is ' , precis
  else
    write ( *, '(a)' ) '  logarithmic encoding (radix=1) has precision'
    write ( *, '(a)' ) '  characterized solely by U1.'
  end if

1650    continue

  if ( fp1 <= ulppls * fp9 * fp9 * t8 ) then
    sdefct = sdefct + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect: '
    write ( *, '(a)' ) '  Precision less than 5 significant decimals,'
    write ( *, '(a)' ) '  which is usually inadequate.'
  end if

  return
end
subroutine random ( x, y, x9, r9 )

!*****************************************************************************80
!
!! RANDOM computes a "somewhat" randomized number.
!
!  Modified:
!
!    19 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Output, real X, the value of (X9+R9)**5.
!
!    Output, real Y, the fractional part of X.
!
!    Input/output, real X9, on output, a "somewhat" randomized number.
!
!    Input, real R9, a number used in the randomization.
!
  implicit none

  real r9
  real x
  real x9
  real y

  x = x9 + r9
  y = x * x
  y = y * y
  x = x * y
  y = x - real ( int ( x ) )

  x9 = y + x * 0.000005E+00

  return
end
subroutine round ( a1, divgrd, fails, flaws, fp1, fp2, fp3, fp9, half, &
  miles, minus1, mulgrd, onemin, r1, r2, r3, radix, sticky, subgrd, &
  ulpmin, ulppls )

!*****************************************************************************80
!
!! ROUND tests rounding in multiplication, division, addition and subtraction.
!
!  Modified:
!
!    18 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    ?, integer MILES, the number of milestone reached so far in testing.
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  real a
  real a1
  real b1
  real b2
  real b9
  real divgrd
  integer fails
  integer flaws
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  real fp3
  real fp9
  real half
  integer miles
  real minus1
  real mulgrd
  real onemin
  real q
  real r1
  real r2
  real r3
  real radix
  real s1
  real sticky
  real subgrd
  real t
  real t5
  real ulpmin
  real ulppls
  real x
  real y
  real y1
  real y2
  real z

  b2 = radix / fp2
  t5 = fp1 + half
  b9 = ( ( radix - fp1) - ulppls ) + fp1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROUND:'
  write ( *, '(a)' ) '  Checking for rounding in multiply, divide, '
  write ( *, '(a)' ) '  add and subtract:'

  r1 = minus1
  r2 = minus1
  r3 = minus1
!
!  Is the radix a power of 2 or 10?
!
  a1 = fp2

  do

    a = radix

    do

      x = a
      a = a / a1

      if ( aint ( a ) /= a ) then
        exit
      end if

    end do

    if ( x == fp1 ) then
      exit
    end if
!
!  Radix is a power of A1; if radix=1 then  a1=2.
!
    if ( fp2 < a1 ) then
      a1 = radix
      exit
    end if

    a1 = fp9 + fp1

  end do
!
!  Is radix a power of 10?
!
!  Unless B is a power of a1 and a1 = 2 or 10.
!
  a = fp1 / a1
  x = a1
  y = a

  do

    z = x * y - half

    if ( z /= half ) then
      fails = fails + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *,2361) x, y, x, x
2361    format(' 1/',e13.5,' = ',e13.5,', and  ', &
               e13.5,'*(1/',e13.5,') differs from  1.')
    end if

    if ( x == radix ) then
      exit
    end if

    x = radix
    y = fp1 / x

  end do

  y2 = fp1 + ulppls
  y1 = fp1 - ulppls
  x = t5 - ulppls
  y = t5 + ulppls
  z = ( x - ulppls ) * y2
  t = y * y1
  z = z - x
  t = t - x
  x = x * y2
  y = ( y + ulppls ) * y1
  x = x - t5
  y = y - t5

  if ( x == fp0 .and. y == fp0 .and. z == fp0 .and. t <= fp0 ) then

    x = ( t5 + ulppls) * y2
    y = t5 - ulppls - ulppls
    z = t5 + ulppls + ulppls
    t = ( t5 - ulppls ) * y1
    x = x - ( z + ulppls )
    sticky = y * y1
    s1 = z * y2
    t = t - y
    y = ( ulppls - y ) + sticky
    z = s1 - ( z + ulppls + ulppls )
    sticky = ( y2 + ulppls ) * y1
    y1 = y2 * y1
    sticky = sticky - y2
    y1 = y1 - half

    if ( x == fp0 .and. y      == fp0 .and. z  == fp0 .and. &
         t == fp0 .and. sticky == fp0 .and. y1 == half ) then
      r1 = fp1
    end if

    if ( x + ulppls        == fp0 .and. y  < fp0 .and. &
         z + ulppls        == fp0 .and. t  < fp0 .and. &
         sticky + ulppls == fp0 .and. y1 < half ) then
      r1 = fp0
    end if

    if ( r1 == fp0 ) then
      write ( *, '(a)' ) '  Multiplication appears to be chopped.'
    else if (r1 == fp1) then
      write ( *, '(a)' ) '  Multiplication appears to be correctly rounded.'
    end if

    if ( r1 - mulgrd == fp1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a)' ) '  Multiplication test is inconsistent;'
      write ( *, '(a)' ) '  Please notify Karpinski!'
    end if

  end if

  if ( r1 == minus1 ) then
    write ( *, '(a)' ) &
      '  Multiplication is neither chopped nor correctly rounded.'
  end if

  miles = 45

  y2 = fp1+ulppls
  y1 = fp1-ulppls
  z = t5+ulppls+ulppls
  x = z/y2
  t = t5-ulppls-ulppls
  y = (t-ulppls)/y1
  z = (z+ulppls)/y2
  x = x-t5
  y = y-t
  t = t/y1
  z = z-(t5+ulppls)
  t = (ulppls-t5)+t

  if ( x <= fp0 .and. y <= fp0 .and. z <= fp0 .and. t <= fp0 ) then

    x = t5/y2
    y = t5-ulppls
    z = t5+ulppls
    x = x-y
    t = t5/y1
    y = y/y1
    t = t-(z+ulppls)
    y = y-z
    z = z/y2
    y1 = (y2+ulppls)/y2
    z = z-t5
    y2 = y1-y2
    y1 = (onemin-ulpmin)/onemin

    if ( x == fp0 .and. y  == fp0 .and. z == fp0   .and. &
         t == fp0 .and. y2 == fp0 .and. y1-half == onemin-half ) then
      r2=fp1
    end if

    if ( x < fp0 .and. y  < fp0 .and. z       < fp0   .and. &
      t < fp0 .and. y2 < fp0 .and. y1-half < onemin-half ) then
      r2=fp0
    end if

    if ( r2 == fp0 ) then
      write ( *, '(a)' ) '  Division appears to be chopped.'
    end if

    if ( r2 == fp1 ) then
      write ( *, '(a)' ) '  Division appears to be correctly rounded.'
    end if

    if ( r2-divgrd == fp1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a)' ) '  Division test is inconsistent;'
      write ( *, '(a)' ) '  Please notify Karpinski!'
    end if

  end if

  if ( r2 == minus1 ) then
    write ( *, '(a)' ) &
    '  Division is neither chopped nor correctly rounded.'
  end if

  b1 = fp1 / radix

  if ( b1 * radix - half /= half) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  radix * (1 / radix)  differs from  1.'
  end if

  miles = 50

  if ( (onemin + ulpmin) - half /= half         .or. &
    (b9 + ulppls) - fp1  /= radix - fp1 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Incomplete carry-propagation in addition.'
  end if

  x = fp1 - ulpmin * ulpmin
  y = fp1 + ulppls * (fp1 - ulppls)
  z = onemin - half
  x = (x - half) - z
  y = y - fp1

  if ( x == fp0 .and. y == fp0 ) then
    r3 = fp0
    write ( *, '(a)' ) '  Add/subtract appears to be chopped.'
  end if

  if ( subgrd /= fp0 ) then

    x = (half + ulppls) * ulppls
    y = (half - ulppls) * ulppls
    x = fp1 + x
    y = fp1 + y
    x = (fp1 + ulppls) - x
    y = fp1 - y

    if ( x == fp0 .and. y == fp0 ) then

      x = (half + ulppls) * ulpmin
      y = (half - ulppls) * ulpmin
      x = fp1 - x
      y = fp1 - y
      x = onemin - x
      y = fp1 - y

      if ( x == fp0 .and. y == fp0 ) then

        r3 = minus1 - fp2 * r3
        write ( *, '(a)' ) '  Add/subtract appears to be correctly rounded.'
  
        if ( r3 - subgrd == fp1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Failure:'
          write ( *, '(a)' ) '  Add/subtract test is inconsistent;'
          write ( *, '(a)' ) '  Please notify Karpinski!'
        end if

      end if

    end if

  end if

  if ( r3 == minus1 ) then
    write ( *, '(a)' ) '  Add/subtract neither chopped nor correctly rounded.'
  end if

  s1 = fp1
  x = fp1+half*(fp1+half)
  y = (fp1+ulppls)*half
  z = x-y
  t = y-x
  sticky = z+t

  if ( sticky /= fp0 ) then
    s1 = fp0
    flaws = flaws + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Flaw:'
    write ( *, '(a,e16.8)' ) '  Nonzero (X-Y) + (Y-X) = ', sticky
    write ( *,2760) x, y
2760    format('  when x = ',e16.8,'  and  y = ',e16.8)
  end if

  sticky = fp0

  if ( mulgrd * divgrd * subgrd < fp1 .or. r1 < fp1 .or. r2 < fp1 .or. &
    r3 < fp1 .or. aint(b2) /= b2 ) go to 2890

        write ( *, '(a)' ) '  Checking for sticky bit:'

  x = (half+ulpmin)*ulppls
        y = half*ulppls
        z = fp1+y
        t = fp1+x
        if ( fp0 < z-fp1 .or. t-fp1 < ulppls) go to 2890
    z = t+y
        y = z-x
        if (z-t < ulppls .or. y-t /= fp0) go to 2890
    x = (half+ulpmin)*ulpmin
        y = half*ulpmin
        z = fp1-y
        t = fp1-x

  if (z-fp1 /= fp0 .or. t-onemin /= fp0) go to 2890

  z=(half-ulpmin)*ulpmin
  t=onemin-z
  q=onemin-y

  if (t-onemin /= fp0 .or. (onemin-ulpmin)-q /= fp0) go to 2890

  z = ( fp1 + ulppls ) * t5
  t = ( t5 + ulppls ) - z + ulppls
  x = fp1 + half / radix
  y = fp1 + radix * ulppls
  z = x * y

  if ( t /= fp0 .or. (x + radix * ulppls) - z /= fp0 ) go to 2890

  if ( radix /= fp2 ) then
    x = fp2 + ulppls
    y = x / fp2
    if ( y - fp1 /= fp0 ) then
      go to 2890
    end if
  end if

  sticky = s1

  if ( sticky == fp1 ) then
    write ( *, '(a)' ) '  Sticky bit appears to be used correctly.'
  end if

2890    continue

  if ( sticky == fp0 ) then
    write ( *, '(a)' ) '  Sticky bit used incorrectly or not at all.'
  end if

  if ( mulgrd * divgrd * subgrd == fp0 .or. &
            r1 < fp0                   .or. &
            r2 < fp0                   .or. &
            r3 < fp0                        ) then
    flaws = flaws + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Flaw:'
    write ( *, '(a)' ) '  Lack of guard digits or failure' 
    write ( *, '(a)' ) '  to correctly round or chop (noted above) '
    write ( *, '(a)' ) '  count as one flaw in the final tally below.'
  end if

  return
end
subroutine small_int ( defect, fails, fp1, fp2, from, ieee, miles, &
  sdefct )

!*****************************************************************************80
!
!! SMALL_INT carries out tests on small integers.
!
!  Modified:
!
!    20 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    ?, integer FROM, the number of the milestone to return to on restart.
!
!    ?, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    ?, integer MILES, the number of milestone reached so far in testing.
!
  implicit none

  real, parameter :: temp12 = 3.0E+00 * 4.0E+00
  real, parameter :: temp20 = 4.0E+00 * ( 4.0E+00 + 1.0E+00 )

  integer defect
  real, parameter :: eight = 4.0E+00 + 4.0E+00
  integer fails
  real, parameter :: five = 4.0E+00 + 1.0E+00
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  integer from
  real, parameter :: half = 1.0E+00 / 2.0E+00
  integer ieee
  integer miles
  real, parameter :: minone = -1.0E+00
  real, parameter :: nine = 3.0E+00 * 3.0E+00
  integer partu
  real radix_temp
  logical restrt
  integer sdefct
  real temp
  real, parameter :: temp240 = temp20 * temp12
  real temp27
  real temp32
  real temp48
  real temp60
  real temp80
  real tempz
  real ulppls_temp
  character ( len = 8 ), parameter :: charz = 'z'

  if (from == 7) go to 951

  if ( from /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SMALL_INT - Error!'
    write ( *, '(a)' ) '  Unrecognized restart milestone.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SMALL_INT'
  write ( *, '(a)' ) '  Small integer tests.'
!
!  Look for some obvious mistakes.
!
  if ( 0.0E+00 + 0.0E+00 /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 0+0=0.'
  end if

  if ( 1.0E+00 - 1.0E+00 /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 1-1=0.'
  end if

  if ( 1.0E+00 <= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 0 < 1.'
  end if

  if ( 1.0E+00 + 1.0E+00 /= 2.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 1+1=2.'
  end if

  temp = 0.0E+00
  tempz = -temp
  if ( tempz == 0.0E+00 ) go to 960
  fails = fails + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Failure:'
  write ( *, '(a)' ) '  Comparison alleges that minus zero, obtained by'
  write ( *, '(a)' ) '  setting X = 0 and then Z = -X, is nonzero!'
!
!  Call to routine to check for partial underflow using minus
!  zero don't really have info on what a unit in the last
!  place is or what the radix is since we haven't gotten to
!  such sophisticated stuff yet, so pick some arbitrary values
!  for now to get us through this next test.
!
  ulppls_temp = 0.001E+00
  radix_temp = 1

951     continue

  restrt = ( from == 7 )

  call partuf ( tempz, charz, defect, fp1, fp2, ieee, miles, partu, &
    radix_temp, restrt, sdefct, ulppls_temp )

960     continue

  if ( 4.0E+00 + 2.0E+00 * (-2.0E+00) /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 4 + 2*(-2) = 0.'
  end if

  if ( ( 4.0E+00 - 3.0E+00 ) - 1.0E+00 /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of (4-3)-1 = 0.'
  end if

  if ( minone + 1.0E+00 /= 0.0E+00 .or. &
    1.0E+00 + minone /= 0.0E+00 .or. &
    minone + abs ( minone ) /= 0.0E+00 .or. &
    minone + minone * minone /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of -1 + 1 = 0.'
  end if

  if ( half + minone + half /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 1/2 - 1 + 1/2 = 0.'
  end if

  miles = 10
  temp27 = nine * 3.0E+00
  temp32 = 4.0E+00 * eight

  if ( temp32 - temp27 - 4.0E+00 - 1.0E+00 /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of   32 - 27 - 4 - 1 = 0.'
  end if

  temp80 = temp240 / 3.0E+00
  temp = temp80 - 4.0E+00 * temp20

  if ( temp /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 240/3 = 80.'
  end if

  temp60 = temp240 / 4.0E+00
  temp = temp60 - five * temp12

  if ( temp /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 240/4 = 60.'
  end if

  temp48 = temp240 / five
  temp = temp48 - 4.0E+00 * temp12

  if ( temp /= 0.0E+00 ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Violation of 240/5 = 48.'
  end if

  if ( fails == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  -1, 0, 1/2 , 1, 2, 3, 4, 5, 9, 27, 32 & 240 are OK.'
  end if

  return
end
subroutine sqrerr ( x, u, j, e5, e7, radix, serous )

!*****************************************************************************80
!
!! SQRERR assesses error in SQRT ( X * X ) - X.
!
!  Modified:
!
!    15 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  real b1
  real e5
  real e6
  real e7
  integer j
  real radix
  logical serous
  real u
  real x

  b1 = 1.0E+00 / radix
  e6 = ( ( sqrt ( x * x ) - x * b1 ) - ( x - x * b1 ) ) / u

  if ( e6 == 0.0E+00 ) then
    return
  end if

  e5 = min ( e5, e6 )
  e7 = max ( e7, e6 )    

  j = j + 1

  if ( .not. serous ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'
    write ( *, '(a,e15.7,a,e15.7,a,e15.7)' ) '  sqrt(', x*x,') - ', x, &
      ' = ', u*e6
    write ( *, '(a)' ) '  instead of correct value of 0.'
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, '(a,e15.7,a,e15.7,a,e15.7)' ) '  sqrt(', x*x,') - ', x, &
      ' = ', u*e6
    write ( *, '(a)' ) '  instead of correct value of 0.'
  end if

  return
end
subroutine sqrtdx ( x, z2, i, d, y2, y, x8, e5, e7, w, radix )

!*****************************************************************************80
!
!! SQRTDX tests the SQRT function.
!
!  Discussion:
!
!    This routine tests if
!
!      sqrt ( d * x ) = sqrt ( (y - 1/2)^2 + x8 / 2 )
!
!    rounds to Y.
!
!  Modified:
!
!    14 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  real d
  real e5
  real e6
  real e7
  real, parameter :: half = 0.5E+00
  integer i
  real radix
  real w
  real x
  real x2
  real x8
  real y
  real y2
  real z2

  if ( x - radix < z2 - radix ) then
    return
  end if

  if ( w - z2 < x - z2 ) then
    return
  end if

  i = i + 1

  x2 =  sqrt ( x * d )
  y2 = ( x2 - z2 ) - ( y - z2 )
  x2 = x8 / ( y - half )
  x2 = x2 - half * x2 * x2

  e6 = ( y2 + half ) + ( half - x2 )
  e5 = min ( e5, e6 )

  e6 = y2 - x2
  e7 = max ( e7, e6 )

  return
end
subroutine square ( a1, defect, fails, fp1, fp2, fp4, fp8, fp9, from, &
  half, ieee, miles, minus1, onemin, precis, r4, radix, sdefct, &
  ulpmin, ulppls, w )

!*****************************************************************************80
!
!! SQUARE tests the square root function.
!
!  Modified:
!
!    19 November 2001
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    ?, integer FROM, the number of the milestone to return to on restart.
!
!    ?, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    ?, integer MILES, the number of milestone reached so far in testing.
!
!    Input, integer RADIX, the computed radix of the machine.
!
  implicit none

  real a1
  real b1
  real b2
  real b9
  integer bad
  real d
  real d4
  integer defect
  real e5
  real e6
  real e7
  integer fails
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  real fp4
  real fp8
  real fp9
  integer from
  real half
  integer i
  integer ieee
  integer j
  integer miles
  real minus1
  integer, parameter :: numtry = 20
  real onemin
  real precis
  real q
  real r4
  real radix
  integer sdefct
  real temp
  real temp1
  real temp2
  real u
  real ulpmin
  real ulppls
  real w
  real x
  real x1
  real x8
  real y
  real y1
  real y2
  real z
  real z1
  real z2
!
!  Trap integer overflows.
!  Other exceptions are trapped by default.
!
  b2 = radix / fp2
  b9 = ( ( radix - fp1 ) - ulppls ) + fp1
  b1 = fp1 / radix

  if ( from == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SQUARE:'
    write ( *, '(a)' ) '  Running tests of square root.'

    miles = 79
    call check_write ( ieee, miles )
    x = fp0
    i = 0

    do

      y = sqrt ( x )

      if ( y /= x .or. y - half /= x - half ) then
        fails = fails + 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Failure:'
        write ( *, '(a,e9.1,a,e15.7)' ) '  SQRT(', x, '), miscalculated as ', y
      end if

      x = -x
      i = i + 1

      if ( i /= 1 ) then
        exit
      end if

    end do

  else if ( from == 79 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  SQRT(-0.0) stops the machine.'
    fails = fails + 1
    i = 2

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SQUARE - Error!'
    write ( *, '(a)' ) '  Unrecognized restart milestone.'
    stop

  end if

  x = fp1
  i = i + 1

  y = sqrt(x)

  if ( y /= x .or. y - half /= x - half ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a,e9.1,a,e15.7)' ) '  SQRT(', x, '), miscalculated as ', y
  end if

  x = -x
  i = i + 1

  x = fp1
  i = i + 1
!
!  Record min and max errors.
!
  e5 = fp0
  e7 = fp0
!
!  Test whether sqrt(x*x) = x.
!
  j = 0

  x = radix
  u = ulppls
  call sqrerr ( x, u, j, e5, e7, radix, .true. )

  x = b1
  u = b1 * ulpmin
  call sqrerr ( x, u, j, e5, e7, radix, .true. )

  x = w
  u = fp1
  call sqrerr ( x, u, j, e5, e7, radix, .true. )

  x = ulpmin
  u = ulpmin * ulpmin
  call sqrerr ( x, u, j, e5, e7, radix, .true. )
!
!  If sqrt has serious defects, then pause.
!
  if ( j /= 0 ) then
    sdefct = sdefct + j
    call check_write ( ieee, miles )
    miles = miles + 1
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Testing if SQRT(X*X) = X for various integers X.'

  j = 0
  x = fp2
  y = radix

  if ( radix /= fp1 ) then

    do

      x = y
      y = radix * x

      if ( numtry <= y - x ) then
        exit
      end if

    end do

  end if

  u = x * ulppls

  bad = 0

  do i = 1, numtry

    x = x + fp1
    call sqrerr ( x, u, j, e5, e7, radix, .false.)

    if ( 0 < j ) then
      bad = bad + j
      exit
    end if

  end do

  defect = defect + bad

  if ( bad == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Found no discrepancies.'
  end if
!
!  Test for monotonicity.
!
  i = -1
  x = b9
  y = radix
  z = radix + radix * ulppls

  bad = 0

  do

    i = i + 1
    x = sqrt ( x )
    q = sqrt ( y )
    z = sqrt ( z )

    if ( q < x .or. z < q ) then
      bad = bad + 1
      defect = defect + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Defect:'
      write ( *, '(a,e15.7)' ) '  SQRT(X) is not monotonic for X near ', y
      exit
    end if

    q = aint ( q + half )

    if ( i <= 0 .and. q * q /= radix ) then
      exit
    end if

    if ( i <= 0 ) then
      y = q
      x = y - ulppls
      z = y + ulppls
    else if ( i == 1 ) then
      y = y * b1
      x = y - ulpmin
      z = y + ulpmin
    else
      exit
    end if

  end do

  if ( bad == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  SQRT has passed a test for monotonicity.'
  end if

3390    continue

  miles = 80
!
!  Test SQRT for accuracy.
!
!  e5 = min{error + 1/2}
!  e7 = max{error - 1/2}
!
  e5 = e5 + half
  e7 = e7 - half
  y = ( sqrt ( fp1 + ulppls ) - fp1 ) / ulppls
  e6 = ( y - fp1 ) + ulppls / fp8
  e7 = max ( e7, e6 )
  e6 = y + ulppls / fp8
  e5 = min ( e5, e6 )
  y = (( sqrt ( onemin ) - ulppls ) - ( fp1 - ulppls ) ) / ulpmin
  e6 = y + ulpmin / fp8
  e7 = max ( e7, e6 )
  e6 = ( y + fp1 ) + ulpmin / fp8
  e5 = min ( e5, e6 )
  i = 0
  u = ulppls
  x = u
!
!  loop
!
  do

    i = i + 1
    y =  sqrt ( ( x + ulpmin + x ) + onemin )
    y = ( ( y - ulppls ) - ((fp1 - ulppls) + x ) ) / u
    z = ( ( ulpmin - x ) + onemin ) * half * x * x / u
    e6 = ( y + half ) + z
    e5 = min ( e5, e6 )
    e6 = ( y - half ) + z
    e7 = max ( e7, e6 )

    if ( i == 4 ) then
      exit
    end if

    if ( i == 2 ) then
      u = ulpmin
      x = -u
    else
      x = u *  sign ( fp1, x ) * aint ( fp8 / ( fp9 *  sqrt ( u ) ) )
    end if

  end do

  miles = 85
  r4 = minus1

  if ( radix == fp1 ) then
    go to 3900
  end if

  write ( *, '(a)' ) '  Testing whether SQRT is rounded or chopped:'

  d = aint ( half + radix ** ( fp1 + precis - aint ( precis ) ) )
!
!   =  b^(1 + fract)  if  p  =  integer  +  fract.
!
  x = d / radix
  y = d / a1

  if ( x /= aint ( x ) .or. y /= aint ( y ) ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a,e15.7)' ) '  Anomalous arithmetic with integers < B^P = ', w
    write ( *, '(a)' ) '  foils testing whether SQRT rounds or chops.'
    go to 3940
  end if

  x = fp0
  z2 = x
  y = fp1
  y2 = y
  z1 = radix - fp1
  d4 = fp4 * d
!
!  loop: for y = 1, 3, 5, ...  maximize y2 = y*y mod 4d .
!
  do

    if ( z2 < y2 ) then

      q = radix
      y1 = y
!
!  If old < new y2, check that gcd(y,b) = 1.
!
      do

        temp = half - q / y1
        temp1 = aint ( temp )

        if ( temp < temp1 ) then
          temp1 = temp1 - fp1
        end if

        x1 = abs ( q + temp1 * y1 )
        q = y1
        y1 = x1

        if ( x1 <= fp0 ) then
          exit
        end if

      end do
!
!  If 1 < gcd(y,b), then skip over y;  else
!
      if ( q <= fp1 ) then
        z2 = y2
        z = y
      end if

    end if
!
!  and gcd(z, radix) = 1
!
    y = y + fp2
    x = x + fp8
    y2 = y2 + x

    if ( .not. (y2 < d4) ) then
      y2 = y2 - d4
    end if
!
!  =  y*y mod 4d
!
    if ( d <= y ) then
      exit
    end if

  end do
!
!  else  0 < z < d  &  z2 = z^2 mod 4d  is maximal.
!
  x8 = d4 - z2
  q = (x8 + z * z) / d4
  x8 = x8 / fp8

  if ( q == aint ( q ) ) then

    do

      x = z1 * z
      x = x - aint(x / radix) * radix

      if ( x == fp1 ) then
        go to 3800
      end if
!
!  with 1 = z*z1 mod b
!
      z1 = z1 - fp1

      if ( z1 <= fp0 ) then
        exit
      end if

    end do

  end if

  fails = fails + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Failure:'
  write ( *, '(a,e15.7)' ) '  Anomalous arithmetic with integers < B^P = ', w
  write ( *, '(a)' ) '  foils testing whether SQRT rounds or chops.'
  go to 3940
!
!   - b/2 <= z1 == 1/z mod b <= b/2
!
3800    continue

  if ( b2 < z1 ) then
    z1 = z1 - radix
  end if
!
!  Loop until d =  b^(p - 1) .
!
  do

    call newd (x, z1, q, z, d, fp1, fp2, half, radix )

    if ( onemin <= ulppls * d ) then
      exit
    end if

  end do

  if ( d * radix - d /= w - d ) then
    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a,e15.7)' ) '  Anomalous arithmetic with integers < B^P = ', w
    write ( *, '(a)' ) '  foils testing whether SQRT rounds or chops.'
    go to 3940
  end if

  z2 = d
  i = 0
!
!  Count how many tests of  sqrt(d*x) = y yield results.
!
    y = d + (fp1 + z) * half
    x = d + z + q
    call sqrtdx (x, z2, i, d, y2, y, x8, e5, e7, w, radix )

    y = d + (fp1 - z) * half + d
    x = d - z + d
    x = x + q + x
    call sqrtdx (x, z2, i, d, y2, y, x8, e5, e7, w, radix )

    call newd (x, z1, q, z, d, fp1, fp2, half, radix )

    if ( d - z2 /= w - z2 ) then
      fails = fails + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a,e15.7)' ) &
        '  Anomalous arithmetic with integers < B^P = ', w
      write ( *, '(a)' ) '  foils testing whether SQRT rounds or chops.'
      go to 3940
    end if

    y = ( d - z2 ) + ( z2 + ( fp1 - z ) * half )
    x = ( d - z2 ) + ( z2 - z + q )
    call sqrtdx (x, z2, i, d, y2, y, x8, e5, e7, w, radix )
    y = (fp1 + z) * half
    x = q
    call sqrtdx (x, z2, i, d, y2, y, x8, e5, e7, w, radix )

    if ( i == 0 ) then
      fails = fails + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a,e15.7)' ) &
        '  Anomalous arithmetic with integers < B^P = ', w
      write ( *, '(a)' ) '  foils testing whether SQRT rounds or chops.'
      go to 3940
    end if

3900    continue

  if ( 0.0E+00 <= e5 .and. e7 <= 0.0E+00 ) then
    r4 = fp1
    write ( *, '(a)' ) '  SQRT appears to be correctly rounded.'
    return
  end if

  if ( e7 + ulppls <= ulppls - half .and. e5 <= half .and. &
    half < e5 + radix ) then
    r4 = fp0
    write ( *, '(a)' ) '  SQRT appears to be chopped.'
    return
  end if

3940    continue

  write ( *, '(a)' ) '  SQRT is neither chopped nor correctly rounded.'
  temp = e5 - half
  temp2 = half + e7
  write ( *,3942) temp, temp2
3942    format('  Observed errors run from  ',e15.7,'  to  ',e15.7, &
       ' ulps.')

  if ( radix * radix <= e7 - e5 ) then
    sdefct = sdefct + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, '(a)' ) '  SQRT gets too many last digits wrong.'
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 )  time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
subroutine underflow ( c1, fp1, fp2, fp3, fp8, from, h1, half, &
  ieee, miles, minpos, nulps, phony0, uflthr )

!*****************************************************************************80
!
!! UNDERFLOW tests underflow threshholds.
!
!  Modified:
!
!    29 March 2003
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, integer FROM, the number of the milestone to return to on restart.
!
!    ?, integer IEEE, a flag whether gradual underflows are doubly rounded.
!
!    ?, integer MILES, the number of milestone reached so far in testing.
!
!    ?, integer NUMTRY, the number of times to try random trials.
!
!  Local Parameters:
!
!    accur  ... flag to indicate success/failure of accuracy tests
!    error  ... count of errors detected testing powers.
!    i      ... scratch for enumerating cases
!    iq     ... temporary for holding integer exponents
!    partu  ... flag to indicate the detecion of partial underflow
!    c      ... 1/(radix^large_integer)
!    epsp1  ... epsilon + 1 (1 + (small integer)* 1 ulp of 1+...)
!    exp2   ... value of e ^ 2
!    h      ... min (1/radix, 1/2)
!    mindif ... minimum positive number found by addition/subtr.
!
  implicit none

  integer accur
  real c
  real c1
  character ( len = 8 ), parameter :: chare0 = ' minpos'
  character ( len = 8 ), parameter :: charz0 = ' phony0'
  real d
  real epsp1
  integer error
  real exp2
  real, parameter :: fp0 = 0.0E+00
  real fp1
  real fp2
  real fp3
  real, parameter :: fp32 = 32.0E+00
  real fp8
  integer from
  real h
  real h1
  real half
  integer i
  integer ieee
  integer iq
  integer miles
  real mindif
  real minpos
  real nulps
  integer, parameter :: numtry = 20
  integer partu
  real phony0
  real q
  real r
  logical restrt
  real t0
  real temp
  real uflthr
  real v9
  real x
  real y
  real y1
  real y2
  real z
  real z9

  common /global/ fails, sdefct, defect, flaws, radix, ulppls, &
    ulpmin, precis, w, mulgrd, divgrd, subgrd, a1, onemin
  integer fails, sdefct, defect, flaws
  real radix, ulppls, &
    ulpmin, precis, w, mulgrd, divgrd, subgrd, a1, onemin

  if ( from /= 0 ) then
!
!  We must be doing a restart.  Figure out where, and go do it
!  must read the log file back in.
!
    read(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
    read(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
    rewind 3
    if (from == 105) go to 4390
    if (from == 106) go to 4410
    if (from == 107) go to 4450
    if (from == 108) phony0 = 0
    if (from == 108) go to 4522
    if (from == 109) go to 4860
    if (from == 115) go to 4631
    if (from == 120) go to 4890
    if (from == 121) go to 4941
    if (from == 122) go to 5011
    if (from == 123) go to 5160
    if (from == 124) go to 5190
    if (from == 125) go to 5175

    if (from == 131) then
      flaws = flaws + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Flaw:'
      write ( *, '(a)' ) '  Underflow trap from exponentiation **.'
      go to 5310
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'UNDERFLOW - Fatal error!'
    write ( *, '(a,i6)' ) '  Unrecognized restart milestone FROM = ', from
    stop

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNDERFLOW:'
  write ( *, '(a)' ) '  Seeking underflow threshold'
  write ( *, '(a)' ) '  and minimum positive number:'
  miles = 105
  call check_write ( ieee, miles )
  d = ulpmin

  if ( precis /= aint ( precis ) ) then

    d = fp1 / radix
    x = precis

    do
      d = d / radix
      x = x - fp1
      if ( x <= fp0 ) then
        exit
      end if
    end do

  end if
!
!  If non-integral precision, we now have D = 1 right shifted by PRECIS
!  digits in base RADIX.
!
!  If integral precision, ULPMIN is this number - pre-computed.
!
  y = fp1
  z = d
!
!  D = a power of  1/RADIX < 1
!
  do
    
    c = y
    y = z

    write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
    write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
    rewind 3

    call check_write ( ieee, miles )
    z = y * y

    if ( y <= z .or. z + z <= z ) then
      exit
    end if

  end do

4390    continue

  miles = 106
  call check_write ( ieee, miles )
  y = c
  z = y * d

  do

    c = y
    y = z
    write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
    write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
    rewind 3
    call check_write ( ieee, miles )
    z = y * d

    if ( y <= z .or. z+z <= z ) then
      exit
    end if

  end do

4410    continue

  miles = 107
  call check_write ( ieee, miles )
  h1 = radix
  if ( h1 < fp2 ) then
    h1 = fp2
  end if

  h = fp1 / h1
!
!  1/h1 = h = min{ 1/radix, 1/2 }
!
  c1 = fp1 / c
  minpos = c
  z = minpos * h
!
!  c = 1/radix^(big integer) << 1 << c1 = 1/c
!
  do

    y = minpos
    minpos = z
    write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
    write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
    rewind 3
    call check_write ( ieee, miles )
    z = minpos * h

    if ( minpos <= z .or. z+z <= z) then
      exit
    end if

  end do

4450    continue

  miles = 108
  call check_write ( ieee, miles )
  uflthr = minpos
  mindif = fp0
  q = fp0
  nulps = ulppls
  epsp1 = fp1 + nulps
  d = c * epsp1

  if ( d <= c ) then

    nulps = radix * ulppls
    epsp1 = fp1 + nulps
    d = c * epsp1
!
!  Multiplication is too crude.
!
    if ( d <= c ) then
      fails = fails + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Failure:'
      write ( *, '(a)' ) '  Multiplication gets too many last digits wrong.'
      t0 = minpos
      y1 = fp0
      phony0 = z
      call check_write ( ieee, miles )
      miles = miles + 1
      go to 4570
    end if

  end if

  t0 = d
  phony0 = t0 * h
  uflthr = fp0

  do

    y1 = t0
    t0 = phony0

    if ( mindif + mindif <= mindif ) then

      y2 = t0 * h1
      mindif = abs ( y1 - y2 )
      q = y1
      if ( uflthr == fp0 .and. y1 /= y2 ) then
        uflthr = y1
      end if

      write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
      write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
      rewind 3
      miles = 108
      call check_write ( ieee, miles )

    end if

    phony0 = t0 * h

    if ( t0 <= phony0 .or. phony0 + phony0 <= phony0 ) then
      exit
    end if

  end do

4522 continue

  miles = 109
  call check_write ( ieee, miles )
!
!  Now 1 >> c=1/radix^(integer) >=  y > minpos=y*h >~ z:=minpos*h >~ 0 ,
!  and 1 >> d=(1+nulps)*c >= uflthr >= q >= y1 > t0:=y1*h >~ phony0:=t0*h >~ 0 ,
!  and uflthr = d/radix^integer is first to violate (uflthr*h)/h=uflthr, 
!  else uflthr=0 ;
!  and q:=uflthr/radix^integer is first with  mindif := |(q*h)/h - q| > 0, 
!  else q=y1.
!
4570    continue

  if ( phony0 == fp0 ) go to 4860
!
!  Test  phony0  for 'phoney-zero' violating 
!    phony0<t0 
!  or
!    phony0<phony0+phony0  ...
!
  write ( *, '(a)' ) ' '
  z = phony0

  if ( phony0 <= fp0 ) then

    fails = fails + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Failure:'
    write ( *, '(a)' ) '  Positive expressions can underflow to an'
    write ( *, '(a)' ) '  allegedly negative value Z0 that '
    write ( *, '(a,e16.8)' ) '  prints out as ', phony0
    x = -phony0

    if ( x <= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Moreover, -Z0, which should then be positive,'
      write ( *, '(a,e16.8)' ) '  is not.  It prints out as ', x
    end if

  else

    flaws = flaws + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Flaw:'
    write ( *, '(a)' ) '  Underflow can stick at an allegedly positive' 
    write ( *, '(a,e16.8)' ) '  value Z0 that prints as ', phony0

  end if

  miles = 115

  write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
  write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
  rewind 3

4631    continue

  restrt = ( from == 115 )

  call partuf ( z, charz0, defect, fp1, fp2, ieee, miles, partu, &
    radix, restrt, sdefct, ulppls )
!
!  End of test for 'phoney-zero'.
!
4860    continue

  miles=120
  write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
  write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
  rewind 3
  call check_write ( ieee, miles )
!
!  as happens on most machines.
!
    if ( c1 * y1 < c1 * y ) then
      epsp1 = h * epsp1
      minpos = t0
    end if
!
!       = least positive number on hp 3000
!
4890    continue

  if ( mindif /= 0 .and. mindif /= minpos ) then

    if ( minpos <= mindif ) then
      defect = defect + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Defect:'
      write ( *, '(a)' ) '  Differences underflow at a higher threshold '
      write ( *, '(a)' ) '  than products.'
    else
      defect = defect + 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Defect:'
      write ( *, '(a)' ) '  Products underflow at a higher threshold'
      write ( *, '(a)' ) '  than differences.'
      if ( phony0 == fp0 ) then
        minpos = mindif
      end if
    end if

  end if
!
!  but not if pseudo-zeros exist.
!
  write ( *, '(a,1pe16.8)' ) '  Smallest positive number found is ', minpos
  z = minpos
  miles = 121
  write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
  write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
  rewind 3

4941    continue

  restrt = ( from == 121 )

  call partuf ( z, chare0, defect, fp1, fp2, ieee, miles, partu, &
    radix, restrt, sdefct, ulppls )

  if ( partu == 1 ) then
    t0 = y
  else
    t0 = minpos
  end if
!
!  for cdc 7600
!
  if ( mindif == fp0 ) then
    i = 3
  else
    i = 4
  end if
!
!  i=1 if mindif=0=uflthr,   
!  i=2 if mindif>0=uflthr,
!
4960    continue

  if ( uflthr == fp0 ) then
    i = i - 2
  end if
!
!  i=3 if mindif=0<uflthr, 
!  i=4 if mindif>0 & uflthr>0
!
  go to (4980, 5090, 5010, 5130),i
!
!       ... case statement
!
4980    uflthr=t0
        if ( c1*q == (c1*y)*epsp1 ) go to 5010
        fails = fails + 1
        uflthr=y
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Failure:'
        write ( *, '(a)' ) '  Either accuracy deteriorates as numbers approach'
        write ( *, '(a)' ) '  a threshold of'
        write ( *,4996)uflthr,c
4996    format(e16.8,' coming down from  ',e16.8,',')
        write ( *, '(a)' ) '  or else multiplication gets too many last '
        write ( *, '(a)' ) '  digits wrong.'
        call check_write ( ieee, miles )
        miles = miles + 1
!
!  Test for  x-z = 0  although  x  /=  z
!
5010    miles = 122
        write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
        write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
        rewind 3
        call check_write ( ieee, miles )
        r =  sqrt ( t0 / uflthr )
        go to 5012
5011    r = fp1
5012    continue

        if ( r <= h ) then
          z = r * uflthr
          x = z * (fp1+r*h*(fp1+h))
        else
          z = uflthr
          x = z * (fp1+h*h*(fp1+h))
        end if

        miles = 123
        write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
        write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
        rewind 3
        call check_write ( ieee, miles )

        if ( x == z .or. x-z /= fp0 ) go to 5160
        flaws = flaws + 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Flaw:'
        write ( *,5055)x,z
5055    format('  x =',e16.8,' is unequal to  z =',e16.8,' ,')
        z9 = x - z
        write ( *,5057) z9
5057    format(' yet  x-z  yields ', e15.7)
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Should this not signal underflow, this is '
        write ( *, '(a)' ) '  a serious defect that causes confusion when'
        write ( *, '(a)' ) '  innocent statements like'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    if ( x == z ) then'
        write ( *, '(a)' ) '      ...'
        write ( *, '(a)' ) '    else'
        write ( *, '(a)' ) '      fp = ( f(x) - f(z) ) / ( x - z )'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  encounter division by zero, although'
        write ( *, '(a)' ) ' '
        write ( *, '(a,e16.8)' ) '    X / Z = 1 + ', (x/z-half)-half
        go to 5160
!
!  End of test for  x-z = 0  &  x  /=  z
!
5090    continue
!
!  case i=2
!  uflthr = 0 < mindif
!
  fails = fails + 1
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Failure:'
  write ( *, '(a)' ) '  Underflow confuses comparison, which alleges'
  write ( *, '(a)' ) '  that Q = Y while denying that |Q-Y| = 0.'
  write ( *, '(a)' ) '  These values print out as:'
  temp = abs ( q - y2 )
  write ( *,5106)q,y2,temp
5106    format(' q =',e16.8,',  y =',e16.8,',  |q-y| =',e16.8,' ,')
  temp = q / y2 - half
  write ( *,5110) temp - half
5110    format(' and  q/y = 1 + ',e16.8)
  uflthr = q
  go to 5010
!
!  case i=4 ;  uflthr > 0  and mindif > 0
!
5130    continue

    if (.not. (q == uflthr .and. mindif == minpos .and. &
          abs ( uflthr - mindif / nulps ) <= mindif ) ) go to 5010

  write ( *, '(a)' ) '  Underflow is gradual; it incurs absolute error = '
  write ( *, '(a)' ) '  (roundoff in underflow threshold) < MINPOS.'
  y = minpos*c1
  y = y * ( 1.5E+00 + ulppls )
  x = c1 * ( fp1 + ulppls )
  y = y / x

  if ( y == minpos ) then
    ieee = 1
  else
    ieee = 0
  end if
!
!  IEEE=1 unless gradual underflows are doubly rounded.
!
5160    continue

  write ( *, '(a)' ) ' '
  write ( *, '(a,e16.8)' ) '  The underflow threshhold is ', uflthr
  write ( *, '(a)' ) '  below which, calculation may suffer larger relative' 
  write ( *, '(a)' ) '  error than merely roundoff.'

  miles = 124

  write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
  write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
  rewind 3

  call check_write ( ieee, miles )

  y2 = ulpmin * ulpmin
  y = y2 * y2
  miles = 125

  write(3) accur, c, epsp1, error, exp2, h, i, iq, mindif, partu
  write(3) d, q, r, t0, v9, x, y, y1, y2, z, z9
  rewind 3

  call check_write ( ieee, miles )

  y2 = y * ulpmin

5175    continue

  if ( uflthr < y2 ) then
    go to 5220
  end if

  if ( minpos < y ) then
    defect = defect + 1
    i = 5
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'
    write ( *, '(a)' ) '  Range is too narrow.'
    write ( *, '(a,i1)' ) '  Underflow for ULPMIN ^', i
    go to 5220
  end if

5190    continue

  sdefct = sdefct + 1
  i=4
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Serious defect:'
  write ( *, '(a)' ) '  Range is too narrow.'
  write ( *, '(a,i1)' ) '  Underflow for ULPMIN ^', i

5220    continue

  miles = 130
  call check_write ( ieee, miles )
  miles = miles + 1

  y = -aint ( half - 240.0E+00 * log ( uflthr ) / log ( h1 ) ) / 240
  y2 = y + y
  write ( *,5240)h1,y
5240    format(' since underflow occurs below the threshold  ='/10x,'(', &
           1pe16.8,')^(',1pe16.8,') ,')
        write ( *,5245)h1,y2
5245    format(' only underflow should afflict the expression'/10x,'(', &
      1pe16.8,')^(',1pe16.8,') ;')

  write ( *, '(a)' ) '  Actually calculating it yields:'
  miles = 131
  call check_write ( ieee, miles )
  v9 = h1 ** (y2)
  write ( *, '(e16.8)' ) v9

  if ( v9 < fp0 .or. (radix+radix*nulps)*uflthr < v9 ) then
    sdefct = sdefct + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Serious defect:'
    write ( *, '(a,e16.8)' ) &
      '  This is not between 0 and underflow threshold  =', uflthr
  else if ( v9 <= uflthr * ( fp1 + nulps ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This computed value is OK.'
  else
    defect = defect + 1
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Defect:'
    write ( *, '(a,e16.8)' ) &
      '  This is not between 0 and underflow threshold  =', uflthr
  end if

5310    continue

  miles = 140
!
!  Calculate exp2 = exp(2) = 7.389056099...
!
  x = fp0
  i = 2
  y = fp2 * fp3
  q = fp0
  accur = 0

  do

    z = x
    i = i + 1
    y = y / ( i + i )
    r = y + q
    x = z + r
    q = ( z - x ) + r

    if ( x <= z ) then
      exit
    end if

  end do

  z = ( 1.5E+00 + fp1 / fp8 ) + x / ( 1.5E+00 * fp32 )
  x = z*z
  exp2 = x * x
  x = onemin
  y = x - ulpmin

  write ( *,5360) exp2
5360    format('  Testing  x^((x+1)/(x-1)) vs. exp(2) = ',e16.8,'  as  x -> 1.')

  do

    do i = 1, numtry

      z = x - ( 1.0E+00 / radix )
      z = ( x + fp1 ) / ( z - ( fp1 - ( 1.0E+00 / radix ) ) )
      q = x**z - exp2

      if ( 240.0E+00 * ulppls < abs ( q ) ) then
        accur = 1
        defect = defect + 1
        temp = + ( x - ( 1.0E+00 / radix ) ) - ( fp1 - ( 1 / radix ) )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Defect:'
        write ( *,5425)temp,z
5425    format('  calculated  (1 + (',e16.8,'))^(',e16.8,')')
        write ( *,5427)q
5427    format('         differs from correct value by  ',e16.8)
        write ( *, '(a)' ) '  This much error may spoil financial' 
        write ( *, '(a)' ) '  calculations involving tiny interest rates.'
        go to 5450
      end if

      z = ( y - x ) * fp2 + y
      x = y
      y = z
      z = fp1 + ( x - onemin ) * ( x - onemin )

      if (z <= fp1 ) then
        exit
      end if

    end do

    if ( fp1 < x ) then
      if ( accur == 0 ) then
        write ( *, '(a)' ) '  Accuracy seems adequate.'
      end if
      exit
    end if

    x = fp1 + ulppls
    y = ulppls + ulppls + x

  end do

5450    continue

  miles=150

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Testing powers Z^Q at four nearly extreme values:'

  error = 0
  z = a1
  iq = int ( half - log ( c ) / log ( a1 ) )

  do

    x = c1
    call cmpxy ( defect, x, y, z, iq, error )
    iq = -iq
    x = c
    call cmpxy ( defect, x, y, z, iq, error )

    if ( z < fp1 ) then
      exit
    end if

    z = 1.0E+00 / a1

  end do

  if ( 0 < error ) then
    write ( *, '(a,i4,a)' ) '  Similar discrepancies have occurred ', error, &
      ' times.'
  end if

  if ( error == 0 ) then
    write ( *, '(a)' ) '  No discrepancies found.'
  else
    call check_write ( ieee, miles )
    miles = miles + 1
  end if

  miles = 160

  return
end
subroutine zeros ( fp0, from, ieee, miles )

!*****************************************************************************80
!
!! ZEROS investigates division by zero.
!
!  Modified:
!
!    17 January 2003
!
!  Author:
!
!    W M Kahan
!
!  Parameters:
!
!    Input, real FP0, the value of 0.
!
!    Input/output, integer FROM, the number of the milestone to return 
!    to on restart.
!
!    Input, integer IEEE, flag whether gradual underflows are doubly rounded.
!
!    Output, integer MILES, the number of milestone reached so far in testing.
!
  implicit none

  real fp0
  real, parameter :: fp1 = 1.0E+00
  integer from
  integer ieee
  integer miles
  real q9

  if ( from == 0 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ZEROS:'
    write ( *, '(a)' ) '  Now determine the error messages and values'
    write ( *, '(a)' ) '  produced when dividing by zero.'

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  About to compute 1/0:'

    miles = 211
    call check_write ( ieee, miles )

    q9 = fp1 / fp0

    write ( *, '(a)' ) ' '
    write ( *, '(a,1pe15.7)' ) '  1/0 = ', q9

    from = 211

  end if

  if ( from == 211 ) then

    miles = 212

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  About to compute 0/0:'

    call check_write ( ieee, miles )
    q9 = fp0 / fp0

    write ( *, '(a)' ) ' '
    write ( *, '(a,1pe15.7)' ) '  0/0 = ', q9

    from = 212

  end if

  if ( from == 212 ) then

    miles = 220

  end if

  return
end

