program main

!*****************************************************************************80
!
!! PCPRB6 demonstrates the jacobian approximation options.
!
!  Modified:
!
!    10 November 1999
!
!  Author:
!
!    John Burkardt
!
!    RWORK(18) is used to determine the size of the finite difference
!    increments used to approximate the jacobian.  The continuation
!    program will set a default value for this, but the user may
!    override it.
!
!    IWORK(1) can be used to check the user supplied jacobian routine
!    against a forward or central difference approximations, and to
!    specify whether the maximum value of the difference, or the
!    entire matrix of differences is to be printed.
!
!    The problem solved is the Freudenstein-Roth function, although
!    the solution procedure is only carried out for five steps.  In
!    this example, we're more interested in the output of the
!    jacobian checker than in the solution of the problem.
!
!    On the second try of the problem, the Jacobian is "polluted".
!    The jacobian checker points this out.  Note, however, that the
!    program is able to compute points on the curve, even with a bad
!    Jacobian.  This is typical.  The good news is that small errors
!    in the jacobian aren't fatal.  The bad news is that you lose
!    quadratic convergence in Newton's method, and may not realize
!    that something's wrong, since you still get answers!
!
  integer, parameter :: nvar = 3
  integer, parameter :: liw = nvar + 29
  integer, parameter :: lrw = 29 + ( 6 + nvar ) * nvar

  external dge_slv
  external dfbad
  external dfroth
  external fxroth
  external pitcon

  double precision fpar(1)
  integer i
  integer ierror
  integer ipar(1)
  integer isave
  integer itry
  integer iwork(liw)
  integer j
  character ( len = 12 ) name
  double precision rwork(lrw)
  double precision xr(nvar)

  call timestamp ( )
  write ( *, * ) ' '
  write ( *, * ) 'PITCON7_PRB6:'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  PITCON test problem'
  write ( *, * ) '  Freudenstein-Roth function'
  write ( *, * ) ' '
  write ( *, '(a,i8)' ) '  Number of equations is ', nvar - 1
  write ( *, '(a,i8)' ) '  Number of variables is ', nvar
  write ( *, * ) ' '
  write ( *, * ) 'This test demonstrates the use of IWORK(1)'
  write ( *, * ) 'and RWORK(18) to approximate the jacobian,'
  write ( *, * ) 'compare the user jacobian to an approximation,'
  write ( *, * ) 'choose forward or centered differences,'
  write ( *, * ) 'choose the size of the difference increment,'
  write ( *, * ) 'print user, approximate jacobian or difference,'
  write ( *, * ) 'print full matrix, or maximum entry.'
  write ( *, * ) ' '

  itry=0
10    continue
  itry=itry+1

  write ( *, * ) ' '

  if(itry==1)then
    write ( *, * ) 'Run 1: standard run for comparison.'
    write ( *, * ) 'Don''t call jacobian checker.'
  elseif(itry==2)then
    write ( *, * ) 'Run 2: run with bad jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value IWORK(1)=-1.'
  elseif(itry==3)then
    write ( *, * ) 'Run 3: run with bad jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-2.'
  elseif(itry==4)then
    write ( *, * ) 'Run 4: run with bad jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-3.'
  elseif(itry==5)then
    write ( *, * ) 'Run 5: run with bad jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-4.'
  elseif(itry==6)then
    write ( *, * ) 'Run 6: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-4.'
    write ( *, * ) 'Use default value of rwork(18)'
  elseif(itry==7)then
    write ( *, * ) 'Run 7: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-4.'
    write ( *, * ) 'Use finite difference increment rwork(18)=0.1'
  elseif(itry==8)then
    write ( *, * ) 'Run 8: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-4.'
    write ( *, * ) 'Use finite difference increment rwork(18)=0.01'
  elseif(itry==9)then
    write ( *, * ) 'Run 9: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-4.'
    write ( *, * ) 'Finite difference increment rwork(18)=0.0001'
  elseif(itry==10)then
    write ( *, * ) 'Run 10: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-5.'
  elseif(itry==11)then
    write ( *, * ) 'Run 11: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-6.'
  elseif(itry==12)then
    write ( *, * ) 'Run 12: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-7.'
  elseif(itry==13)then
    write ( *, * ) 'Run 13: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-8.'
  elseif(itry==14)then
    write ( *, * ) 'Run 14: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-9.'
  elseif(itry==15)then
    write ( *, * ) 'Run 15: run with good jacobian.'
    write ( *, * ) 'Call jacobian checker at third step.'
    write ( *, * ) 'Use check value iwork(1)=-10.'
  end if

  write ( *, * ) ' '

  iwork(1:liw) = 0
  rwork(1:lrw) = 0.0
!
!  IWORK(1)=0 ; This is a startup
!  IWORK(2)=2 ; Use index 2 for first parameter
!  IWORK(3)=0 ; Program may choose index
!  IWORK(4)=0 ; Update jacobian every newton step
!  IWORK(5)=3 ; Seek target values for index 3
!  IWORK(6)=1 ; Seek limit points in index 1
!  IWORK(7)=0 ; small amount of output
!  IWORK(9)=0 ; Use user's jacobian routine
!
  iwork(1)=0
  iwork(2)=2
  iwork(3)=0
  iwork(4)=0
  iwork(5)=3
  iwork(6)=1
  iwork(7)=0
  iwork(9)=0
!
!  RWORK(1)=0.00001; Absolute error tolerance
!  RWORK(2)=0.0001 ; Relative error tolerance
!  RWORK(3)=0.01   ; Minimum stepsize
!  RWORK(4)=20.0   ; Maximum stepsize
!  RWORK(5)=0.3    ; Starting stepsize
!  RWORK(6)=1.0    ; Starting direction
!  RWORK(7)=1.0    ; Target value (Seek solution with X(3)=1)
!
  rwork(1)=0.00001
  rwork(2)=0.0001
  rwork(3)=0.01
  rwork(4)=20.0
  rwork(5)=0.3
  rwork(6)=1.0
  rwork(7)=1.0
!
!  Set the parameter that determines the size of the increment used for
!  finite difference approximations.
!
  if(itry==7)then
    rwork(18)=0.1
  elseif(itry==8)then
    rwork(18)=0.01
  elseif(itry==9)then
    rwork(18)=0.0001
  else
    rwork(18)=0.0
  end if
!
!  Set starting point.
!
  xr(1)=15.0
  xr(2)=-2.0
  xr(3)=0.0

  write ( *, * ) ' '
  write ( *, * ) 'Step  Type of point     X(1)         X(2)         X(3)'
  write ( *, * ) ' '
  i=0
  name='Start point  '
  write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)
!
!  Take another step.  This loop normally runs for 40 steps, but we'll cut it
!  to just 5.
!
  do 40 i=1,5
!
!  Check jacobian on third step
!
    if(i==3.and.itry>1)then

      isave=iwork(1)

      if(itry==2)then
        iwork(1)=-1
      elseif(itry==3)then
        iwork(1)=-2
      elseif(itry==4)then
        iwork(1)=-3
      elseif(itry==5)then
        iwork(1)=-4
      elseif(itry==10)then
        iwork(1)=-5
      elseif(itry==11)then
        iwork(1)=-6
      elseif(itry==12)then
        iwork(1)=-7
      elseif(itry==13)then
        iwork(1)=-8
      elseif(itry==14)then
        iwork(1)=-9
      elseif(itry==15)then
        iwork(1)=-10
      else
        iwork(1)=-4
      end if

      write ( *, * ) ' '
      write ( *, * ) 'Turning on jacobian check option!'
      write ( *, * ) ' '
    end if

    if ( 2 <= itry .and. itry <= 5 ) then
      call pitcon(dfbad,fpar,fxroth,ierror,ipar,iwork,liw, &
        nvar,rwork,lrw,xr,dge_slv)
    else
      call pitcon(dfroth,fpar,fxroth,ierror,ipar,iwork,liw, &
        nvar,rwork,lrw,xr,dge_slv)
    end if

    if(iwork(1)==1)then
      name='Corrected   '
    elseif(iwork(1)==2)then
      name='Continuation'
    elseif(iwork(1)==3)then
      name='Target point'
    elseif(iwork(1)==4)then
      name='Limit point '
    elseif(iwork(1)<0)then
      name='Jacobian'
    end if
!
!  After Jacobian check, restore value of IWORK(1).
!
    if(i==3.and.itry>1)then
      iwork(1)=isave
    end if

    write(*,'(1x,i3,2x,a12,2x,3g14.6)')i,name,(xr(j),j=1,nvar)

    if(iwork(1)==3)then
      write ( *, * ) ' '
      write ( *, * ) 'We have reached the point we wanted.'
      go to 50
    end if

    if ( ierror /= 0 ) then
      write ( *, * ) ' '
      write ( *, * ) 'PITCON returned an error code:'
      write ( *, * ) 'IERROR = ', ierror
      write ( *, * ) ' '
      write ( *, * ) 'The computation failed.'
      go to 50
    end if

40      continue

50    continue

  if(itry<15)go to 10
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PITCON7_PRB6:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine fxroth ( nvar, fpar, ipar, x, f )
!
!*****************************************************************************80
!
!! FXROTH evaluates the function F(X) at X.
!
!
!  The function has the form:
!
!    ( X1 - ((X2-5.0)*X2+2.0)*X2 - 13.0 + 34.0*(X3-1.0)  )
!    ( X1 + ((X2+1.0)*X2-14.0)*X2 - 29.0 + 10.0*(X3-1.0) )
!
  integer nvar
!
  double precision f(*)
  double precision fpar(*)
  integer ipar(*)
  double precision x(nvar)
!
  f(1) = x(1) - ( ( x(2) - 5.0 ) * x(2) + 2.0 ) * x(2)- 13.0 &
         + 34.0 * ( x(3) - 1.0 )

  f(2) = x(1) + ( ( x(2) + 1.0 ) * x(2) - 14.0 ) * x(2) - 29.0 &
         + 10.0 * ( x(3) - 1.0 )

  return
end
subroutine dfroth ( nvar, fpar, ipar, x, fjac )

!*****************************************************************************80
!
!! DFROTH evaluates the Jacobian J(X) at X.
!
!  The jacobian has the form:
!
!    ( 1.0   (-3.0*X(2)+10.0)*X(2)- 2.0   34.0  )
!    ( 1.0   ( 3.0*X(2)+ 2.0)*X(2)-14.0   10.0  )
!
  integer nvar

  double precision fjac(nvar,nvar)
  double precision fpar(*)
  integer ipar(*)
  double precision x(nvar)

  fjac(1,1) = 1.0
  fjac(1,2) = ( - 3.0 * x(2) + 10.0 ) * x(2) - 2.0
  fjac(1,3) = 34.0

  fjac(2,1) = 1.0
  fjac(2,2) = ( 3.0 * x(2) + 2.0 ) * x(2) - 14.0
  fjac(2,3) = 10.0

  return
end
subroutine dfbad ( nvar, fpar, ipar, x, fjac )

!*****************************************************************************80
!
!! DFBAD computes an incorrect value of the the jacobian.
!
!  Discussion:
!
!    This routine is simply used to demonstrate how errors in the
!    jacobian can affect the calculation, and can be caught by using
!    the jacobian checking option.
!
  integer nvar

  double precision fjac(nvar,nvar)
  double precision fpar(1)
  integer ipar(1)
  double precision x(nvar)
!
!  Get the correct values.
!
  call dfroth ( nvar, fpar, ipar, x, fjac )
!
!  Now perturb them.
!
  fjac(1,1) = fjac(1,1) + 0.125
  fjac(1,2) = fjac(1,2) + 0.250
  fjac(2,2) = fjac(2,2) + 0.500

  return
end
