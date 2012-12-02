program main

!*****************************************************************************80
!
!! MAIN is the main program for CYSTAL_QED.
!
!  Discussion:
!
!    CRYSTAL_QED seeks parameters that minimize a crystallization cost functional.
!
!  Discussion:
!
!    The DQED package is used to solve the bounded and constrained least
!    squares problem.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none
!
!  Parameters that don't depend on anything else.
!
  integer ( kind = 4 ), parameter :: liopt = 17
  integer ( kind = 4 ), parameter :: lropt = 1
  integer ( kind = 4 ), parameter :: mcon = 6
  integer ( kind = 4 ), parameter :: mequa = 1
  integer ( kind = 4 ), parameter :: mvars = 7
  integer ( kind = 4 ), parameter :: nt = 5
!
!  First order parameters.
!
  integer ( kind = 4 ), parameter :: ldfj = mcon+mequa
  integer ( kind = 4 ), parameter :: nall = mcon+2*mvars+nt+1
!
!  Second order parameters.
!
  integer ( kind = 4 ), parameter :: liwork = 3*mcon+9*mvars+4*nt+nall+11
  integer ( kind = 4 ), parameter :: lwork = nall*nall+4*nall+mvars*nt+33*mvars &
    +mequa*nt+mvars*nt+14*nt+9*mcon+26+3*nall+2

  real ( kind = 8 ) bl(mvars+mcon)
  real ( kind = 8 ) bu(mvars+mcon)
  real ( kind = 8 ) dnrm2
  external dqed_evaluate
  real ( kind = 8 ) fj(ldfj,mvars+1)
  real ( kind = 8 ) fnorm
  external func
  real ( kind = 8 ) fx(mcon+mequa)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) ind(mvars+mcon)
  integer ( kind = 4 ) iopt(liopt)
  integer ( kind = 4 ) iopti
  integer ( kind = 4 ) ivars
  integer ( kind = 4 ) iwork(liwork)
  integer ( kind = 4 ) nvars
  real ( kind = 8 ) ropt(lropt)
  real tarray(2)
  real ( kind = 8 ) temp
  real time1
  real time2
  real ( kind = 8 ) value
  real ( kind = 8 ) work(lwork)
  real ( kind = 8 ) xpar(mvars)

  call timestamp ( )
!
!  Read the system CPU clock at the starting time.
!
  call cpu_time ( time1 )
!
!  Say hello.
!
  call hello ( liwork, lwork )
!
!  Tell DQED the size of the work arrays.
!
  iwork(1) = lwork
  iwork(2) = liwork
!
!  Set the initial parameter values.
!
  ivars = 0
!
!  These are Hui Zhang's original crucible shape parameters.
!
  if ( ivars == 0 ) then

    nvars = 7
    xpar(1) = 0.015D+00
    xpar(2) = 0.03D+00
    xpar(3) = 0.046D+00
    xpar(4) = 0.07D+00
    xpar(5) = 0.10D+00
    xpar(6) = 0.14D+00
    xpar(7) = 0.18D+00
!
!  These crucible shape parameters represent the optimum solution
!  when the nodes are monotonically constrained.
!
  else if ( ivars == 1 ) then

    nvars = 7
    xpar(1) = 0.064D+00
    xpar(2) = 0.064D+00
    xpar(3) = 0.064D+00
    xpar(4) = 0.25D+00
    xpar(5) = 0.25D+00
    xpar(6) = 0.25D+00
    xpar(7) = 0.25D+00
!
!  These crucible shape parameters cause the grid routine to fail,
!  by generating degenerate control volumes.
!
  else if ( ivars == 2 ) then

    nvars = 7
    value = 0.30D+00
    xpar(1) = 3.0D+00*value/5.0D+00
    xpar(2) = value
    xpar(3) = value-(value-0.25D+00)/5.0D+00
    xpar(4) = value-2.0D+00*(value-0.25D+00)/5.0D+00
    xpar(5) = value-3.0D+00*(value-0.25D+00)/5.0D+00
    xpar(6) = value-4.0D+00*(value-0.25D+00)/5.0D+00
    xpar(7) = value-4.5D+00*(value-0.25D+00)/5.0D+00

  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Initial parameter values:'
  write ( *, '(a)' ) ' '
  do i = 1, nvars
    write ( *, '(i6,g14.6)' ) i, xpar(i)
  end do
!
!  Set the constraints on the variables.
!  My first guess is to just restrain all the variables to be
!  between 0 and 0.25.
!
  if ( ivars == 0 .or. ivars == 1 .or. ivars == 2 ) then

    ind(1) = 3
    bl(1) = 0.0D+00
    bu(1) = 0.25D+00

    ind(2) = 3
    bl(2) = 0.0D+00
    bu(2) = 0.25D+00

    ind(3) = 3
    bl(3) = 0.0D+00
    bu(3) = 0.25D+00

    ind(4) = 3
    bl(4) = 0.0D+00
    bu(4) = 0.25D+00

    ind(5) = 3
    bl(5) = 0.0D+00
    bu(5) = 0.25D+00

    ind(6) = 3
    bl(6) = 0.0D+00
    bu(6) = 0.25D+00

    ind(7) = 3
    bl(7) = 0.0D+00
    bu(7) = 0.25D+00
!
!  The following 6 constraints force monotonicity.
!
     ind(8) = 1
     bl(8) = 0.0D+00
     bu(8) = 0.0D+00

     ind(9) = 1
     bl(9) = 0.0D+00
     bu(9) = 0.0D+00

     ind(10) = 1
     bl(10) = 0.0D+00
     bu(10) = 0.0D+00

     ind(11) = 1
     bl(11) = 0.0D+00
     bu(11) = 0.0D+00

     ind(12) = 1
     bl(12) = 0.0D+00
     bu(12) = 0.0D+00

     ind(13) = 1
     bl(13) = 0.0D+00
     bu(13) = 0.0D+00

   end if
!
!  Set IOPTI, the switch which chooses
!    just one function call (for checking),
!    just one jacobian call (for checking),
!    or optimization.
!
  iopti = 0

  write ( *, '(a)' ) ' '
  if ( iopti == 0 ) then
    write ( *, '(a)' ) 'Just call the function evaluator once.'
  else if ( iopti == 1 ) then
    write ( *, '(a)' ) 'Just call the jacobian approximator once.'
  else
    write ( *, '(a)' ) 'Optimize the problem.'
  end if
!
!  If IOPTI = 0, just call the function evaluator once.
!
  if ( iopti == 0 ) then

    call func(fx,iopt,mcon,mequa,nvars,ropt,xpar)

    fnorm = dnrm2 ( mequa, fx(mcon+1), 1 )
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) 'Two-norm of residual is ', fnorm
!
!  If IOPTI = 1, just call the jacobian approximator once.
!
  else if ( iopti == 1 ) then

    call diffor(fj,func,fx,iopt,ldfj,mcon,mequa,nvars,ropt,xpar)
!
!  If IOPTI = 2, then carry out an optimization.
!
  else

    igo = 0
!
!  IOPT(1) = 2 means change the value of ITMAX to IOPT(2).
!
    iopt(1) = 2
    iopt(2) = 10
    write ( *, '(a)' ) ' '
    write ( *, * ) 'Maximum number of minimization steps is ',iopt(2)
!
!  IOPT(3) = 99 means no more options.
!
    iopt(3) = 99

    call dqed ( dqed_evaluate, mequa,nvars,mcon,ind,bl,bu,xpar,fj,ldfj, &
      fnorm,igo,iopt,ropt,iwork,work)

    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'DQED return flag IGO = ', igo
    write ( *, '(a,i6)' ) 'DQED real work array requirement  = ', iwork(1)
    write ( *, '(a,i6)' ) 'DQED integer work array requirement  = ', iwork(2)
    write ( *, '(a,g14.6)' ) 'Two-norm of residual is ', fnorm
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DQED returned the parameter values:'
    write ( *, '(a)' ) ' '
    do i = 1,nvars
      write ( *, '(i6,g14.6)' ) i,xpar(i)
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) 'with the function value ', fnorm

    call func(fx,iopt,mcon,mequa,nvars,ropt,xpar)

    fnorm = dnrm2(mequa,fx(mcon+1),1)
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I check the output function value!'
    write ( *, '(a,g14.6)' ) 'I get a function value ', fnorm
  end if
!
!  Having completed the desired task, read the system clock and stop.
!
  call cpu_time ( time2 )

  write ( *, '(a)' ) 'CRYSTAL_QED:'
  temp = time2-time1
  write ( *, '(a,g14.6,a)' ) '  Total elapsed CPU time  =  ',temp,' seconds,'
  temp = temp/60.0
  write ( *, '(a,g14.6,a)' ) '                          =  ',temp,' minutes.'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CRYSTAL_QED:'
  write ( *, '(a)' ) ' Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine dqed_evaluate ( xpar, fj, ldfj, igo, iopt, ropt )

!*****************************************************************************80
!
!! DQED_EVALUATE evaluates certain functions being treated by DQED.
!
!  Discussion:
!
!    DQED_EVALUATE also evaluates partial derivatives.
!
!    The user has NVARS variables XPAR(I), and is trying to minimize
!    the square root of the sum of the squares of MEQUA functions
!    F(I)(XPAR), subject to MCON constraints which have the form
!
!      BL(I) < =  G(I)(XPAR) <= BU(I)
!
!    where either the left or right bounding inequality may be dropped.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real ( kind = 8 ) XPAR(*).
!    XPAR is an array of length NVARS, containing the values of the
!    independent variables at which the functions and partial
!    derivatives should be evaluated.
!
!    Output, real ( kind = 8 ) FJ(LDFJ,NVARS+1).
!
!    If IGO is nonzero, then partial derivatives must
!    be placed in the first NVARS columns of FJ, as
!    follows:
!
!      Rows I  =  1 to MCON, and columns J = 1 to NVARS
!      should contain the partials of the I-th constraint
!      function G(I) with respect to the J-th variable.
!
!      Rows I = MCON+1 to MCON+MEQUA, and columns J = 1 to NVARS
!      should contain the partials of the (I-MCON)-th nonlinear
!      function F(I-MCON) with respect to the J-th variable.
!
!    Regardless of the value of IGO, column NVARS+1 must be
!    assigned the values of the constraints and functions,
!    as follows:
!
!      Rows I  =  1 to MCON, column J = NVARS+1, should contain
!      the value of G(I).
!
!      Rows I = MCON+1 to MCON+MEQUA, column J = NVARS+1, should
!      contain the value of F(I-MCON).
!
!    Input, integer ( kind = 4 ) LDFJ.
!    LDFJ is the leading dimension of FJ, which must be at least
!    MCON+MEQUA.
!
!    Input/output, integer ( kind = 4 ) IGO.
!    On input, IGO tells the user whether the partial derivatives
!    are needed.
!
!      0, the partial derivatives are not needed.
!      nonzero, the partial derivatives are needed.
!
!    On output, the user may reset the input value of IGO if one
!    of two situations is encountered:
!
!      99, the functions, constraints, or partial derivatives
!      could not be evaluated at the given input point X.  Request
!      that DQED reject that point, and try a different one.
!
!      Any other value, abort the run.
!
!    Input, integer ( kind = 4 ) IOPT(*), the integer option array.
!
!    Input, real ( kind = 8 ) ROPT(*), the real option array.
!
  implicit none

  integer ( kind = 4 ) ldfj
  integer ( kind = 4 ), parameter :: mcon = 6
  integer ( kind = 4 ), parameter :: mequa = 1
  integer ( kind = 4 ), parameter :: nvars = 7

  real ( kind = 8 ) fj(ldfj,nvars+1)
  external func
  real ( kind = 8 ) fx(mcon+mequa)
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) iopt(*)
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) xpar(nvars)

  if ( igo /= 0 ) then
    call diffor(fj,func,fx,iopt,ldfj,mcon,mequa,nvars,ropt,xpar)
  end if

  call func(fj(1,nvars+1),iopt,mcon,mequa,nvars,ropt,xpar)

  return
end
subroutine dqedhd ( x, fj, ldfj, igo, iopt, ropt )

!*****************************************************************************80
!
!! DQEDHD evaluates functions and derivatives for DQED.
!
!  Discussion:
!
!    For our purposes, this is a dummy routine, and is not called.
!    The compiler is happier when it is here, though.
!
!    The user problem has MCON constraint functions,
!    MEQUA least squares equations, and involves NVARS
!    unknown variables.
!
!    When this subprogram is entered, the general (near)
!    linear constraint partial derivatives, the derivatives
!    for the least squares equations, and the associated
!    function values are placed into the array FJ.
!    all partials and functions are evaluated at the point
!    in X.  then the subprogram returns to the calling
!    program unit. typically one could do the following
!    steps:
!
!    if ( igo /= 0 ) then
!      place the partials of the i-th constraint function with respect to
!      variable j in the array fj(i,j), i = 1,...,mcon, j=1,...,nvars.
!
!    place the values of the i-th constraint equation into fj(i,nvars+1).
!
!    if ( igo /= 0 ) then
!      place the partials of the i-th least squares equation with respect
!      to variable j in the array fj(i,j), i = 1,...,mequa, j = 1,...,nvars.
!
!    place the value of the i-th least squares equation into fj(i,nvars+1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
  implicit none

  integer ( kind = 4 ) ldfj
  integer ( kind = 4 ), parameter :: nvars = 8

  real ( kind = 8 ) fj(ldfj,nvars+1)
  integer ( kind = 4 ) igo
  integer ( kind = 4 ) iopt(*)
  integer ( kind = 4 ) mcon
  integer ( kind = 4 ) mequa
  integer ( kind = 4 ) mode
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) x(nvars)
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQEDHD - Fatal error!'
  write ( *, '(a)' ) '  This is a dummy routine.'
  stop

end
subroutine func ( fx, iopt, mcon, mequa, nvars, ropt, xpar )

!*****************************************************************************80
!
!! FUNC communicates between the optimizer and the constitutive routines.
!
!  Discussion:
!
!    FUNC receives a set of parameter values as input, solves the
!    resulting constitutive equations, and evaluates the derived
!    cost functional, whose value is returned to the calling routine.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) IOPT(*), the integer option array.
!
!    Input, integer ( kind = 4 ) MCON.
!    MCON is the number of constraints imposed on the independent
!    variables that compose a feasible solution.
!
!    Input, integer ( kind = 4 ) MEQUA.
!    MEQUA is the number of components of the nonlinear function.
!
!    Input, integer ( kind = 4 ) NVARS.
!    NVARS is the number of parameters or independent variables
!    upon which the function depends.
!
!    Input, real ( kind = 8 ) ROPT(*).
!    ROPT is the real option array.
!
!    Input, real ( kind = 8 ) XPAR(*).
!    XPAR is an array of length NVARS, containing the values of the
!    parameters.
!
  implicit none
!
  integer ( kind = 4 ), parameter :: maxbot = 20
  integer ( kind = 4 ) mcon
  integer ( kind = 4 ) mequa
  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: nk = 14
  integer ( kind = 4 ), parameter :: ns = 10
  integer ( kind = 4 ) nvars

  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) areal
  real ( kind = 8 ) areas
  real ( kind = 8 ) areat
  real ( kind = 8 ) b
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) b3jbl(ni)
  real ( kind = 8 ) birad
  real ( kind = 8 ) bo
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cfo
  real ( kind = 8 ) cinc
  real ( kind = 8 ) cinco
  real ( kind = 8 ) cmu
  real ( kind = 8 ) cost
  real ( kind = 8 ) cvn
  real ( kind = 8 ) damax
  real ( kind = 8 ) delt
  real ( kind = 8 ) dtm
  real ( kind = 8 ) epsad
  real ( kind = 8 ) epsil
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) fks
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fmax(ns)
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fnu
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) fx(mcon+mequa)
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hf
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icost
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) inturb
  integer ( kind = 4 ) iopt(1)
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipref
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) izone
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) jpref
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lastt
  logical lblk(ns)
  logical lconv
  logical lortho
  logical lsolve(ns)
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nbot
  integer ( kind = 4 ), save :: ncall  =  0
  integer ( kind = 4 ) ndt
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npc
  integer ( kind = 4 ) nsolve(ns)
  integer ( kind = 4 ) ntimes(ns)
  real ( kind = 8 ) orth
  real ( kind = 8 ) pr
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) ra
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) ropt(*)
  real ( kind = 8 ) rpr
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sigt
  real ( kind = 8 ) smax
  real ( kind = 8 ) smooth
  real ( kind = 8 ) ssum
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) tal
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) tas
  real ( kind = 8 ) tend
  real ( kind = 8 ) tf
  real ( kind = 8 ) tinit
  character ( len = 25 ) title(ns)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) tw
  real ( kind = 8 ) vave
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) vort(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xbot(maxbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) xpar(nvars)
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ybot(maxbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ylen
!
  ncall = ncall+1
!
!  Initialize the data.
!
  call inidat(ae1,ae2,ak1,ak2,area,b,b1jbl,b2jbl,b3jbl, &
    birad,bo,cappa,cd,ce1, &
    ce2,cfo,cinc,cmu,cost,cvn,delt,dtm,epsad,epsil,ewall,f,fcsl, &
    fks,fksl,fma,fn,fnu,fo,fr,frsl,gam,gamt,grash,hamag,heta,hf, &
    hksi,icost,icrys,inturb,iplot,ipref,iprint,jcrys, &
    jpref,l0,l1,last,lastt,lblk,lortho,lsolve,m0,m1,maxbot,mode, &
    nbot,ndt,ni,nj,nk,np,npc,ns,nsolve,ntimes,orth, &
    pr,ra,rdtm, &
    re,recb,rect,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk, &
    sigma,sigt,smooth,stel,stes,tal,tanca,tanca2,tas,tend,tf, &
    tinit,title,tnow,tw,vave,vol,vort,x,xbot,xc,xlen,y,ybot,yc,ylen)
!
!  TEMPORARY
!
  do i = 1,7
    xbot(i+1) = xpar(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FUNC:'
  write ( *, '(a)' ) '  Input parameters:'
  do i = 1,7
    write ( *, '(i6,g14.6)' ) i,xpar(i)
  end do
!
!  Print out data.
!
  if ( ncall == 1 ) then
    call prdat(b,birad,bo,cappa,cfo,cvn,delt,dtm,fcsl,fksl,fma, &
      fnu,fr,frsl,grash,icost,icrys,inturb,iprint, &
      jcrys,l0,last,lastt,m0,mode,nbot,ns,nsolve,ntimes,orth, &
      pr,ra,rdtm,recb,rect,rhocon,smooth,stel,stes,tanca,tanca2, &
      tend,tf,tinit,title,tw,vave,xlen,ylen)
  end if
!
!  Generate the initial grid XC, YC.
!
  call inigrd(icrys,jcrys,l0,m0,nbot,xbot,xc,xlen,ybot,yc,ylen)
!
!  Adapt the grid.
!
  do k = 1, 10
    call adapt(cvn,epsad,icrys,iprint,jcrys,l1,m1,nbot, &
      orth,smooth,xbot,xc,ybot,yc)
  end do
!
!  Once XC and YC are determined, compute the control volume areas.
!
  call doarea(area,areal,areas,areat,icrys,jcrys,l0,m0,ni,nj,xc,yc)
!
!  Each time iteration begins at this point.
!
10    continue
!
!  Begin the iterative solution of the state equations in the
!  solid zone.
!
  izone = 1
  l1 = l0
  m1 = jcrys
!
!  Only the temperature (variable 5) needs to be solved for.
!
  lsolve(1) = .false.
  lsolve(2) = .false.
  lsolve(3) = .false.
  lsolve(4) = .false.
  lsolve(5) = .true.
  lsolve(6) = .false.
  lsolve(7) = .false.
  lsolve(8) = .false.
  lsolve(9) = .false.
  lsolve(10) = .false.

  do iter = 0,last

    lconv = .true.
!
!  Set the X and Y coordinates of the primary nodes from XC, YC,
!  the coordinates of the corner nodes.
!
!  Note that this is only done for the left portion of the region!
!
    call setx ( l1, m1, ni, nj, x, xc, y, yc )
!
!  Compute various geometric quantities required for computing
!  derivatives.
!
    call setgeo(ae1,ae2,ak1,ak2,heta,hksi,l1,m1,mode,ni,nj, &
      r,vol,x,xc,y,yc)
!
!  Estimate pressure and momentum at control volume interfaces.
!
    if ( ndt == 0 ) then
      if ( iter == 0 ) then

        call setp(fjeta,fjksi,heta,hksi,l1,m1,ni,nj,f(1,1,3))

        call setru(heta,hksi,l1,m1,ni,nj,rho,rueta,ruksi, &
          f(1,1,1),f(1,1,2),x,y)

      end if
    end if

    call setup(ae1,ae2,ak1,ak2,b1jbl,b2jbl,b3jbl,birad, &
      cappa,cd,ce1,ce2,cmu,epsil,ewall,f, &
      fcsl,fjeta,fjksi,fksl,fma,fmax,fn,fo,frsl,gamt,grash,hamag, &
      heta,hksi,inturb,icrys,iter,izone,jcrys,l1,lblk,lconv,lortho, &
      lsolve,m1,mode,nf,np,npc,nsolve,ntimes,pr,r,rdtm,re,recb, &
      rect,res,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk,sigt, &
      smax,ssum,stel,stes,tal,tas,tf,tw,vol,x,xc,y,yc)

    if ( iprint > 0 ) then
      call output(cfo,iter,izone,res,smax,ssum,f(1,1,5), &
        tnow,f(1,1,1),f(1,1,6))
    end if

    if ( lconv .and. iter >= 5)go to 20

  end do

  if ( iprint > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNC - Warning!'
    write ( *, '(a)' ) '  The solid iteration has not converged'
    write ( *, '(a,i6,a)' ) '  after ',iter,' iterations.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Solid norm and max relative change:'
    write ( *, '(a)' ) ' '
    do i = 1,ns
      if ( lsolve(i) ) then
        fmax(i) = damax ( l0*m0, f(1,1,i), 1 )
        write(*,'(i3,2x,a25,2g14.6)')i,title(i),fmax(i),res(i)
      end if
    end do
  end if
!
!  Zone 2 (LIQUID) calculation
!
20    continue

  do i = 1,icrys
    do j = 1,jcrys
      f(i,j,5) = fo(i,j,5)
      f(i,j,9) = fo(i,j,9)
    end do
  end do

  izone = 2
  l1 = icrys
  m1 = m0

  do i = 1,icrys
    ruksi(i,1) = 0.0D+00
    ruksi(i,m0) = 0.0D+00
    rueta(i,2) = 0.0D+00
    rueta(i,m0) = 0.0D+00
  end do

  do j = 1,m0
    rueta(1,j) = 0.0D+00
    rueta(icrys,j) = 0.0D+00
    ruksi(2,j) = 0.0D+00
    ruksi(icrys,j) = 0.0D+00
  end do

  if ( inturb == 0 ) then
    lsolve(1) = .true.
    lsolve(2) = .true.
    lsolve(3) = .true.
    lsolve(4) = .true.
    lsolve(5) = .true.
    lsolve(6) = .false.
    lsolve(7) = .false.
    lsolve(8) = .false.
    lsolve(9) = .false.
    lsolve(10) = .false.
  else
    lsolve(1) = .true.
    lsolve(2) = .true.
    lsolve(3) = .true.
    lsolve(4) = .true.
    lsolve(5) = .true.
    lsolve(6) = .false.
    lsolve(7) = .true.
    lsolve(8) = .true.
    lsolve(9) = .false.
    lsolve(10) = .false.
  end if

  do iter = 1,last

    lconv = .true.
!
!  Set the turbulent viscosity.
!
    if ( inturb == 0 ) then
      do i = 2,l1-1
        do j = 2,m1-1
          gamt(i,j) = 0.0D+00
        end do
      end do
    else if ( inturb == 1 ) then
      do i = 2,l1-1
        do j = 2,m1-1
          gamt(i,j) = (1.0D+00-relax(12))*gamt(i,j) &
            +relax(12)*cmu*rho(i,j)*f(i,j,7)**2/f(i,j,8)
        end do
      end do
    end if
!
!  Set the X, Y values.
!
    call setx(l1,m1,ni,nj,x,xc,y,yc)
!
!  Compute various geometric quantities required for computing
!  derivatives.
!
    call setgeo(ae1,ae2,ak1,ak2,heta,hksi,l1,m1,mode,ni,nj, &
      r,vol,x,xc,y,yc)
!
!  Set the right hand side of certain flux boundary conditions.
!
    call setup(ae1,ae2,ak1,ak2,b1jbl,b2jbl,b3jbl,birad, &
      cappa,cd,ce1,ce2,cmu,epsil,ewall,f, &
      fcsl,fjeta,fjksi,fksl,fma,fmax,fn,fo,frsl,gamt,grash,hamag, &
      heta,hksi,inturb,icrys,iter,izone,jcrys,l1,lblk,lconv, &
      lortho, &
      lsolve,m1,mode,nf,np,npc,nsolve,ntimes,pr,r,rdtm,re,recb, &
      rect,res,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk,sigt, &
      smax,ssum,stel,stes,tal,tas,tf,tw,vol,x,xc,y,yc)
!
!  Compute the stream function PSI.
!
    f(2,2,10) = 0.0D+00

    do i = 2,icrys

      if ( i > 2 ) then
        f(i,2,10) = f(i-1,2,10)-rueta(i-1,2)*ae1(i-1,2) &
          -0.5D+00*(ruksi(i-1,1)+ruksi(i,1))*ae2(i-1,2)
      end if

      do j = 3,m1-1
        t1 = (rueta(i,j-1)+rueta(i,j))
        t2 = (rueta(i-1,j-1)+rueta(i-1,j))
        t3 = 0.5D+00*(hksi(i,j-1)*t2+hksi(i-1,j-1)*t1)/(hksi(i,j-1)+hksi(i-1,j-1))
        f(i,j,10) = f(i,j-1,10)+ruksi(i,j-1)*ak1(i,j-1)-t3*ak2(i,j-1)
      end do

    end do

    f(1,1,10) = f(2,2,10)

    do j = 2,m0
      f(1,j,10) = f(2,j,10)
    end do

    do i = 2,l0
      f(i,1,10) = f(i,2,10)
    end do
!
!  Optional printout.
!
    if ( iprint > 0 ) then
      call output(cfo,iter,izone,res,smax,ssum,f(1,1,5),tnow,f(1,1,1),f(1,1,6))
    end if

    if ( lconv)go to 30

  end do

  if ( iprint > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FUNC - Warning!'
    write ( *, '(a)' ) '  The liquid iteration has not converged'
    write ( *, '(a,i6,a)' ) '  after ',iter,' iterations.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Liquid norms and max relative change:'
    write ( *, '(a)' ) ' '
    do i = 1,ns
      if ( lsolve(i) ) then
        fmax(i) = damax(l0*m0,f(1,1,i),1)
        write(*,'(i3,2x,a25,2g14.6)')i,title(i),fmax(i),res(i)
      end if
    end do
  end if
!
!  We have computed the state variables for the current time.
!
30    continue
!
!  Compute the vorticity VORT.
!
  call vortic ( heta, hksi, icrys, l0, m0, f(1,1,1), f(1,1,2), vort, xc, yc )
!
!  Evaluate the cost function integrand at the current time.
!
  cinco = cinc
  call setcst(area,cinc,gam,icost,icrys,jcrys,l0,m0,ni,nj,f(1,1,5), &
    tf,f(1,1,1),f(1,1,2),vave,vort,xc,yc)
!
!  Update the total cost, by adding the estimated contribution
!  from the current time interval.
!
!  TEMPORARY
!
  if ( ndt > 0 ) then
!       cost = cost+0.5D+00*delt*(cinco+cinc)
    cost = cost+0.5D+00*dtm*(cinco+cinc)
  end if
!
!  If we've reached the end time, then write out final data,
!  possibly save a restart file, and stop.
!
  if ( ndt >= lastt ) then

    if ( iprint > 0 ) then
      call pmod(ipref,jcrys,jpref,l1,m1,f(1,1,3))
    end if

    cost = cost/(sqrt(areal)*(tend-tinit))
!
!  If requested, write out plot data.
!  Sadly, it is necessary to call SETX to generate the X, Y arrays
!  for the entire region.  Previous calls only generate them for
!  a portion of the region.
!
    if ( iplot == 1 ) then

      call setx(l0,m0,ni,nj,x,xc,y,yc)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEY, ABOUT TO CALL RSWRIT.'

      call rswrit(b1jbl,b2jbl,b3jbl,cost,f(1,1,9),gamt, &
        icrys,jcrys,l0,m0,nbot,f(1,1,3),f(1,1,4),f(1,1,10), &
        rueta,ruksi,f(1,1,5),f(1,1,8),f(1,1,7),tnow,f(1,1,1), &
        f(1,1,2),vort,f(1,1,6),x,xbot,xc,y,ybot,yc)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEY, ABOUT TO CALL WRTEC.'

      call wrtec(gamt,l1,m1,f(1,1,3),f(1,1,10),f(1,1,5),f(1,1,8), &
        f(1,1,7),f(1,1,1),f(1,1,2),f(1,1,6),x,y)

    end if

    write ( *, '(a,g14.6)' ) 'FUNC: COST =        ',cost

    do i = 1,mcon
      fx(i) = xpar(i+1)-xpar(i)
    end do

    fx(mcon+1) = cost

    return

  end if
!
!  Update the number of steps, and the current time.
!
  ndt = ndt+1
  tnow = tinit+ndt*dtm
!
!  Save a copy of the current data.
!
  do i = 1,l0
    do j = 1,m0
      do k = 1,ns
        fo(i,j,k) = f(i,j,k)
      end do
    end do
  end do

  call movgrd(b1jbl,b2jbl,bo,delt,fksl,fr,frsl,icrys, &
    iprint,jcrys,l0,m0,f(1,1,3),pr,re,stel,tanca,tanca2,xc, &
    xlen,yc)

  l1 = l0
  m1 = m0
  call adapt(cvn,epsad,icrys,iprint,jcrys,l1,m1,nbot, &
    orth,smooth,xbot,xc,ybot,yc)
!
!  Once XC and YC are determined, compute the control volume areas.
!
  call doarea(area,areal,areas,areat,icrys,jcrys,l0,m0,ni,nj,xc,yc)

  go to 10
end
subroutine adapt ( cvn, epsad, icrys, iprint, jcrys, l1, m1, nbot, &
  orth, smooth, xbot, xc, ybot, yc )

!*****************************************************************************80
!
!! ADAPT executes MAGG, the multizone adaptive grid generation.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ) nbot
  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  real ( kind = 8 ) a1
  real ( kind = 8 ) a11
  real ( kind = 8 ) a12
  real ( kind = 8 ) a2
  real ( kind = 8 ) a21
  real ( kind = 8 ) a22
  real ( kind = 8 ) a3
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) b3
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) c3
  real ( kind = 8 ) cvn
  real ( kind = 8 ) det
  real ( kind = 8 ) dot
  real ( kind = 8 ) dx
  real ( kind = 8 ) dxc(ni,nj)
  real ( kind = 8 ) dxmax
  real ( kind = 8 ) dy
  real ( kind = 8 ) dyc(ni,nj)
  real ( kind = 8 ) epsad
  real ( kind = 8 ) gn
  real ( kind = 8 ) gt
  real ( kind = 8 ) gx
  real ( kind = 8 ) gy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) nin
  real ( kind = 8 ) orth
  real ( kind = 8 ) pb(ni,nj)
  real ( kind = 8 ) pc1
  real ( kind = 8 ) pc2
  real ( kind = 8 ) pd(ni,nj)
  real ( kind = 8 ) ratio
  real ( kind = 8 ) res1
  real ( kind = 8 ) res2
  real ( kind = 8 ) rmax
  real ( kind = 8 ) s1
  real ( kind = 8 ) s2
  real ( kind = 8 ) scale
  real ( kind = 8 ) smooth
  real ( kind = 8 ) ssmot
  real ( kind = 8 ) ssweg
  real ( kind = 8 ) vmag2
  real ( kind = 8 ) wn
  real ( kind = 8 ) wt
  real ( kind = 8 ) wx
  real ( kind = 8 ) wy
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xdisp
  real ( kind = 8 ) xij
  real ( kind = 8 ) xjac
  real ( kind = 8 ) xn
  real ( kind = 8 ) xnn
  real ( kind = 8 ) xpp
  real ( kind = 8 ) xt
  real ( kind = 8 ) xtn
  real ( kind = 8 ) xtt
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ydisp
  real ( kind = 8 ) yij
  real ( kind = 8 ) yn
  real ( kind = 8 ) ynn
  real ( kind = 8 ) ypp
  real ( kind = 8 ) yt
  real ( kind = 8 ) ytn
  real ( kind = 8 ) ytt
!
  dxmax = 0.0D+00

  do i = 1,l1
    do j = 1,m1
      dxc(i,j) = 0.0D+00
      dyc(i,j) = 0.0D+00
    end do
  end do

  do i = 2,l1
    do j = 2,m1

      if ( i < icrys ) then
        pc1 = 0.95D+00*((i-(2+icrys)/2)/10.0D+00)**2+0.05D+00
      else
        pc1 = 0.8D+00*((i-l1)/12.0)**2+0.2D+00
      end if

      if ( j <= jcrys ) then
        pc2 = 0.9D+00*((j-(2+jcrys)/2)/10.0D+00)**2+0.05D+00
      else
        pc2 = 0.9D+00*((j-(m1+jcrys)/2)/10.0D+00)**2+0.05D+00
      end if

      pb(i,j) = pc1*pc2

      pd(i,j) = pc1*pc2

    end do
  end do
!
!  This DO loop iteration is an iterative solution of a linear system.
!
  do nin = 1,100

    ssmot = 0.0D+00
    ssweg = 0.0D+00
    rmax = 0.0D+00

    if ( iprint > 0 ) then

      if ( nin == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ADAPT: Iter   DXMAX     XC(20,20) YC(20,20)' // &
          ' XC(3,M1)  YC(3,M1)  RMAX'
      end if

      if ( mod(nin,50) == 0 ) then
        write(*,'(6x,i5,6g10.3)') &
          nin,dxmax,xc(20,20),yc(20,20),xc(3,m1),yc(3,m1),rmax
      end if
    end if

    do j = 3,m1-1
      do i = 3,l1-1

        xt = 0.5D+00*(xc(i+1,j)-xc(i-1,j))
        xn = 0.5D+00*(xc(i,j+1)-xc(i,j-1))
        yt = 0.5D+00*(yc(i+1,j)-yc(i-1,j))
        yn = 0.5D+00*(yc(i,j+1)-yc(i,j-1))

        xjac = xt*yn-xn*yt
!
!  Refer to table 2.1, page 14, of Hui Zhang's thesis, for
!  a description of these coefficients.
!
        a1 = -2.0D+00*smooth*(xt*yt+xn*yn)*(xn*xn+yn*yn)/xjac**3
        a2 = 4.0D+00*smooth*(xt*yt+xn*yn)*(xt*xn+yt*yn)/xjac**3
        a3 = -2.0D+00*smooth*(xt*yt+xn*yn)*(xt*xt+yt*yt)/xjac**3

        b1 = 2.0D+00*smooth*(yt*yt+yn*yn)*(xn*xn+yn*yn)/xjac**3
        b2 = -4.0D+00*smooth*(yt*yt+yn*yn)*(xt*xn+yt*yn)/xjac**3
        b3 = 2.0D+00*smooth*(yt*yt+yn*yn)*(xt*xt+yt*yt)/xjac**3

        c1 = 2.0D+00*smooth*(xt*xt+xn*xn)*(xn*xn+yn*yn)/xjac**3
        c2 = -4.0D+00*smooth*(xt*xt+xn*xn)*(xt*xn+yt*yn)/xjac**3
        c3 = 2.0D+00*smooth*(xt*xt+xn*xn)*(xt*xt+yt*yt)/xjac**3

        a1 = a1+orth*pd(i,j)*xn*yn
        a2 = a2+orth*pd(i,j)*(xt*yn+xn*yt)
        a3 = a3+orth*pd(i,j)*xt*yt

        b1 = b1+orth*pd(i,j)*xn*xn
        b2 = b2+2.0D+00*orth*pd(i,j)*(2.0D+00*xt*xn+yt*yn)
        b3 = b3+orth*pd(i,j)*xt*xt

        c1 = c1+orth*pd(i,j)*yn*yn
        c2 = c2+2.0D+00*orth*pd(i,j)*(xt*xn+2.0D+00*yt*yn)
        c3 = c3+orth*pd(i,j)*yt*yt

        a1 = a1-2.0D+00*cvn*pb(i,j)*xn*yn
        a2 = a2+2.0D+00*cvn*pb(i,j)*(xt*yn+xn*yt)
        a3 = a3-2.0D+00*cvn*pb(i,j)*xt*yt

        b1 = b1+2.0D+00*cvn*pb(i,j)*yn*yn
        b2 = b2-4.0D+00*cvn*pb(i,j)*yt*yn
        b3 = b3+2.0D+00*cvn*pb(i,j)*yt*yt

        c1 = c1+2.0D+00*cvn*pb(i,j)*xn*xn
        c2 = c2-4.0D+00*cvn*pb(i,j)*xt*xn
        c3 = c3+2.0D+00*cvn*pb(i,j)*xt*xt
!
!  Now compute terms from the derivatives of the weight functions.
!
        gt = 0.5D+00*(pd(i+1,j)-pd(i-1,j))
        gn = 0.5D+00*(pd(i,j+1)-pd(i,j-1))
        gx = orth*0.5D+00*(gt*yn-gn*yt)*(xt*xn+yt*yn)**2/xjac
        gy = orth*0.5D+00*(gn*xt-gt*xn)*(xt*xn+yt*yn)**2/xjac

        wt = 0.5D+00*(pb(i+1,j)-pb(i-1,j))
        wn = 0.5D+00*(pb(i,j+1)-pb(i,j-1))
        wx = cvn*xjac*(yn*wt-yt*wn)
        wy = cvn*xjac*(xt*wn-xn*wt)
!
!  (RES1, RES2) is the residual that should be driven to zero.
!
        xtt = xc(i+1,j)-2.0D+00*xc(i,j)+xc(i-1,j)
        xnn = xc(i,j+1)-2.0D+00*xc(i,j)+xc(i,j-1)
        ytt = yc(i+1,j)-2.0D+00*yc(i,j)+yc(i-1,j)
        ynn = yc(i,j+1)-2.0D+00*yc(i,j)+yc(i,j-1)
        xtn = 0.25D+00*(xc(i+1,j+1)+xc(i-1,j-1)-xc(i-1,j+1)-xc(i+1,j-1))
        ytn = 0.25D+00*(yc(i+1,j+1)+yc(i-1,j-1)-yc(i-1,j+1)-yc(i+1,j-1))

        res1 = b1*xtt+b2*xtn+b3*xnn+a1*ytt+a2*ytn+a3*ynn+gx+wx
        res2 = a1*xtt+a2*xtn+a3*xnn+c1*ytt+c2*ytn+c3*ynn+gy+wy

        if ( abs(res1)+abs(res2) > rmax ) then
          rmax = abs(res1)+abs(res2)
        end if

        a11 = -2.0D+00*(b1+b3)
        a12 = -2.0D+00*(a1+a3)
        a22 = -2.0D+00*(c1+c3)
        a21 = a12

        det = a11*a22-a12*a21
!
!  (DX, DY) is the pointwise solution of the linear system.
!
        dx = (res2*a21-res1*a22)/det
        dy = (res1*a12-res2*a11)/det
!
!  Now look at how movement of XC(I,J), YC(I,J) affects the neighbors
!
!    (I+1,J-1)     (I+1,J)       (I+1,J+1)
!
!    (I,  J-1)     (I,  J)       (I,  J+1)
!
!    (I-1,J-1)     (I-1,J)       (I-1,J+1)
!
!  to determine the linear factor SCALE that will multiply (DX,DY).
!
        ratio = 0.25D+00
        xij = xc(i,j)
        yij = yc(i,j)

        xdisp = xc(i+1,j)-xij
        ydisp = yc(i+1,j)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i+1,j+1)-xij
        ydisp = yc(i+1,j+1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i,j+1)-xij
        ydisp = yc(i,j+1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i-1,j+1)-xij
        ydisp = yc(i-1,j+1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i-1,j)-xij
        ydisp = yc(i-1,j)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i-1,j-1)-xij
        ydisp = yc(i-1,j-1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i,j-1)-xij
        ydisp = yc(i,j-1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        xdisp = xc(i+1,j-1)-xij
        ydisp = yc(i+1,j-1)-yij
        vmag2 = xdisp*xdisp+ydisp*ydisp
        dot = dx*xdisp+dy*ydisp
        ratio = max(dot/vmag2,ratio)

        scale = min(0.5D+00,0.25D+00/ratio)

        dxc(i,j) = scale*dx
        dyc(i,j) = scale*dy

        ssmot = ssmot+(xt*xt+xn*xn+yt*yt+yn*yn)/xjac
        ssweg = ssweg+cvn*pb(i,j)*xjac**2

      end do
    end do
!
!  Move the points.
!
    do j = 3,m1-1
      do i = 3,l1-1

        if ( i /= icrys ) then
          xc(i,j) = xc(i,j)+dxc(i,j)
        end if

        if ( j /= jcrys ) then
          yc(i,j) = yc(i,j)+dyc(i,j)
        end if

      end do
    end do

    do i = icrys+1,l1-1
      yc(i,jcrys) = 0.4D+00+0.02D+00*sin(xc(i,jcrys)*20.0D+00*3.14159D+00)
    end do
!
!  Handle the nodes along the bottom row, from 10 positions to the right
!  of the crystal, to the right hand wall.
!
    do j = jcrys+10,m1-1

      ratio = 0.25D+00
      call findp(j,xc(3,j),yc(3,j),xpp,ypp,xc,yc)
      dx = xpp-xc(2,j)
      dy = ypp-yc(2,j)

      xdisp = xc(2,j+1)-xc(2,j)
      ydisp = yc(2,j+1)-yc(2,j)
      vmag2 = xdisp*xdisp+ydisp*ydisp
      dot = dx*xdisp+dy*ydisp
      ratio = max(dot/vmag2,ratio)
      s1 = vmag2

      xdisp = xc(2,j-1)-xc(2,j)
      ydisp = yc(2,j-1)-yc(2,j)
      vmag2 = xdisp*xdisp+ydisp*ydisp
      dot = dx*xdisp+dy*ydisp
      ratio = max(dot/vmag2,ratio)
      s2 = max(vmag2,s1)

      scale = min(0.3D+00,0.25D+00/ratio)
      s1 = (scale*dx)**2+(scale*dy)**2

      if ( s1 < s2 ) then
        xc(2,j) = xc(2,j)+scale*dx
        yc(2,j) = yc(2,j)+scale*dy
      end if

      call cubic(nbot,yc(2,j),ybot,xc(2,j),xbot)

      yc(l1,j) = yc(l1-1,j)

    end do
!
!  Handle the nodes along the bottom row, from the left axis
!  of symmetry, to 9 positions beyond the crystal.
!
    do j = 2,jcrys+9

      yc(2,j) = yc(3,j)

      call cubic(nbot,yc(2,j),ybot,xc(2,j),xbot)

      yc(l1,j) = yc(l1-1,j)

    end do

    do i = 3,l1
      xc(i,2) = xc(i,3)
      xc(i,m1) = xc(i,m1-1)
      if ( xc(i,m1-1) <= (xc(2,m1)+0.002D+00* real ( i - 2, kind = 8 ) ) ) then
        xc(i,m1) = xc(2,m1)+0.002D+00* real ( i - 2, kind = 8 )
      end if
    end do
!
!  Compute the maximum movement of all corner nodes.
!
    dxmax = abs(dxc(3,3))
    do j = 3,m1-1
      do i = 3,l1-1
        dxmax = max(dxmax,abs(dxc(i,j)))
        dxmax = max(dxmax,abs(dyc(i,j)))
      end do
    end do
!
!  Copy values to the dummy corner nodes with I = 1 or J=1.
!
    xc(1,1) = xc(2,2)
    yc(1,1) = yc(2,2)

    do j = 2,m1
      xc(1,j) = xc(2,j)
      yc(1,j) = yc(2,j)
    end do

    do i = 2,l1
      xc(i,1) = xc(i,2)
      yc(i,1) = yc(i,2)
    end do
!
!  Check for convergence.
!
    if ( dxmax <= epsad)return

  end do

  return
end
subroutine cubic ( nbot, ynew, ybot, xnew, xbot )

!*****************************************************************************80
!
!! CUBIC constructs a cubic spline through data.
!
!  Discussion:
!
!    CUBIC is given NBOT points (YBOT,XBOT) which lie along the curve
!    defining the bottom of the crucible.
!
!    CUBIC is also given a value YNEW which lies between YBOT(1) and
!    YBOT(NBOT).
!
!    CUBIC constructs a cubic spline through the data, and evaluates it at
!    YNEW, returning the value as XNEW.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) NBOT, the number of points to be used to
!    define the bottom of the crucible.
!
  implicit none
!
  integer ( kind = 4 ) nbot
  integer ( kind = 4 ), parameter :: ni = 64
!
  real ( kind = 8 ) dome
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ist
  real ( kind = 8 ) p2
  real ( kind = 8 ) p3
  real ( kind = 8 ) shift
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tm(ni)
  real ( kind = 8 ) tt(ni)
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xl(ni)
  real ( kind = 8 ) xnew
  real ( kind = 8 ) y1
  real ( kind = 8 ) y2
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yl(ni)
  real ( kind = 8 ) ymin
  real ( kind = 8 ) ynew
!
  if ( ynew < ybot(1) .or. ybot(nbot).lt.ynew ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CUBIC - Fatal error!'
    write ( *, '(a,g14.6)' ) '  YNEW  =  ',ynew
    write ( *, '(a,g14.6)' ) '  outside of range YBOT(1) = ',ybot(1)
    write ( *, '(a,g14.6)' ) '  through YBOT(NBOT) = ',ybot(nbot)
    stop
  end if

  do i = 1,ni
    xl(i) = 0.0D+00
    yl(i) = 0.0D+00
    tt(i) = 0.0D+00
    tm(i) = 0.0D+00
  end do

  yl(1) = ybot(1)+(ybot(1)-ybot(3))
  yl(2) = ybot(2)+(ybot(1)-ybot(3))

  do i = 3,nbot+2
    yl(i) = ybot(i-2)
    xl(i) = xbot(i-2)
  end do

  yl(nbot+3) = yl(nbot+1)+(yl(nbot+2)-yl(nbot))
  yl(nbot+4) = yl(nbot+2)+(yl(nbot+2)-yl(nbot))

  shift = (xl(3)-xl(4))/(yl(3)-yl(4))-(xl(4)-xl(5))/(yl(4)-yl(5))

  xl(2) = xl(3)+(yl(2)-yl(3))*(shift+(xl(3)-xl(4))/(yl(3)-yl(4)))
  xl(1) = xl(2)+(yl(1)-yl(2))*(shift+(xl(2)-xl(3))/(yl(2)-yl(3)))

  shift = (xl(nbot+2)-xl(nbot+1))/(yl(nbot+2)-yl(nbot+1)) &
    -(xl(nbot+1)-xl(nbot))/(yl(nbot+1)-yl(nbot))

  xl(nbot+3) = xl(nbot+2)+(yl(nbot+3)-yl(nbot+2))* &
    (shift+(xl(nbot+2)-xl(nbot+1))/(yl(nbot+2)-yl(nbot+1)))

  xl(nbot+4) = xl(nbot+3)+(yl(nbot+4)-yl(nbot+3))* &
    (shift+(xl(nbot+3)-xl(nbot+2))/(yl(nbot+3)-yl(nbot+2)))

  do i = 1,nbot+3
    tm(i) = (xl(i+1)-xl(i))/(yl(i+1)-yl(i))
  end do

  do i = 3,nbot+2
    dome = abs(tm(i+1)-tm(i))+abs(tm(i-1)+tm(i-2))
    if ( dome < 1.0D-10) then
      tt(i) = 0.0D+00
    else
      tt(i) = (abs(tm(i+1)-tm(i))*tm(i-1)+abs(tm(i-1)-tm(i-2))*tm(i))/dome
    end if
  end do
!
!  Find the node YL(IST) which is nearest to YNEW.
!
  ymin = abs(yl(3)-ynew)
  ist = 3

  do i = 3,nbot+2
    if ( abs(yl(i)-ynew) < ymin ) then
      ist = i
      ymin = abs(yl(i)-ynew)
    end if
  end do

  if ( (yl(ist)-ynew) > 0.0D+00 ) then
    y1 = yl(ist-1)
    x1 = xl(ist-1)
    y2 = yl(ist)
    x2 = xl(ist)
    t1 = tt(ist-1)
    t2 = tt(ist)
  else
    y1 = yl(ist)
    x1 = xl(ist)
    y2 = yl(ist+1)
    x2 = xl(ist+1)
    t1 = tt(ist)
    t2 = tt(ist+1)
  end if
!
!  Evaluate the spline at YNEW.
!
  p2 = (3.0D+00*(x2-x1)/(y2-y1)-2.0D+00*t1-t2)/(y2-y1)
  p3 = (t1+t2-2.0*(x2-x1)/(y2-y1))/(y2-y1)**2

  xnew = x1+t1*(ynew-y1)+p2*(ynew-y1)**2+p3*(ynew-y1)**3

  return
end
subroutine diflow ( acof, diff, flow )
!
!*****************************************************************************80
!
!! DIFLOW computes the convection-diffusion coefficient ACOF.
!
!
!  Discussion:
!
!    DIFLOW is given
!    FLOW, the mass velocity RHO*U, and DIFF, the value of GAMMA/DELX.
!    Several schemes are available, but currently the power law is used.
!
!    See Patankar, chapter 5.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  real ( kind = 8 ) acof
  real ( kind = 8 ) diff
  real ( kind = 8 ) flow

  acof = diff
  if ( diff == 0.0D+00) return
!
!  Power Law.
!
!  Define
!    FE = rho*u
!    DE = gamma/delx
!    Then the coefficient is
!
!    AE  =  DE * Max(0, (1-0.1*Abs(FE)/DE)**5) + Max(0,-FE)
!
  acof = diff*max(1.0D-10,(1.0D+00-0.1D+00*abs(flow/diff))**5)
  if ( flow < 0.0D+00) acof = acof-flow

  return
end
subroutine doarea(area,areal,areas,areat,icrys,jcrys,l0,m0,ni,nj,xc,yc)
!
!*****************************************************************************80
!
!! DOAREA calculates control volume areas.
!
!
!  Discussion:
!
!    DOAREA is given (XC,YC), the locations of the "corners" of the
!    control volumes, and calculates the area of each control volume.
!
!    DOAREA uses the fact that the area of a polygon which is bounded
!    by the N points (X(I),Y(I)) can be computed by:
!
!      AREA  =  0.5 * SUM (I=1 to N) X(I) * (Y(I+1)-Y(I-1))
!
!    where Y(0) should be replaced by Y(N), and Y(N+1) by Y(1).
!
!    Using this formula, AREA may come out negative, depending on
!    whether the nodes are given in clockwise or counterclockwise order.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Output, real ( kind = 8 ) AREA(NI,NJ).
!    AREA contains the area of the control volume (I,J) which is
!    bounded by the corner nodes [I,J], [I+1,J], [I,J+1], and [I+1,J+1].
!
!    Output, real ( kind = 8 ) AREAL, the area of the liquid region.
!
!    Output, real ( kind = 8 ) AREAS.
!    AREAS is the area of the solid or crystal region.
!
!    Output, real ( kind = 8 ) AREAT, the total area.
!
!    Input, integer ( kind = 4 ) ICRYS.
!    ICRYS specifies the end of the crystal in the I array direction,
!    and in the vertical coordinate direction.
!
!    Input, integer ( kind = 4 ) JCRYS.
!    JCRYS specifies the end of the crystal in the J array direction,
!    and in the horizontal coordinate direction.
!
!    Input, integer ( kind = 4 ) L0.
!    L0 is the extent of the grid in the vertical or "I" coordinate.
!
!    Input, integer ( kind = 4 ) M0.
!    M0 is the extent of the grid in the horizontal or "J" coordinate.
!
!    Input, integer ( kind = 4 ) NI.
!    NI is the maximum number of grid points in the "I" or vertical
!    direction.
!
!    Input, integer ( kind = 4 ) NJ.
!    NJ is the maximum number of grid points in the "J" or horizontal
!    direction.
!
!    Input, real ( kind = 8 ) XC(NI,NJ).
!    XC contains the X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Input, real ( kind = 8 ) YC(NI,NJ).
!    YC contains the Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj

  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) areal
  real ( kind = 8 ) areas
  real ( kind = 8 ) areat
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) nbad
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) yc(ni,nj)

  do i = 1,l0
    do j = 1,m0

      if ( i == 1 .or. i == l0.or.j == 1.or.j == m0 ) then
        area(i,j) = 0.0D+00
      else

        area(i,j) = 0.5D+00*( &
           xc(i,  j)*  (yc(i+1,j)  -yc(i,  j+1)) &
          +xc(i+1,j)*  (yc(i+1,j+1)-yc(i,  j)) &
          +xc(i+1,j+1)*(yc(i,  j+1)-yc(i+1,j)) &
          +xc(i,  j+1)*(yc(i,  j)  -yc(i+1,j+1)))

      end if

    end do
  end do
!
!  Check for illegal zero length sides.
!
  nbad = 0
  do i = 2,l0-1
    do j = 2,m0-1
      if ( xc(i,j) == xc(i+1,j) .and. yc(i,j) == yc(i+1,j) ) then
        nbad = nbad+1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETARE - Fatal error!'
        write ( *, '(a,2i6)' ) '  Zero length side for cell I,J:',i,j
        write ( *, '(a,2g14.6)' ) &
          'XC,YC(I,J) = XC,YC(I+1,J)=',xc(i,j),yc(i,j)
      else if ( xc(i+1,j) == xc(i+1,j+1) .and. yc(i+1,j) == yc(i+1,j+1) ) then
        nbad = nbad+1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETARE - Fatal error!'
        write ( *, '(a,2i6)' ) '  Zero length side for cell I,J:',i,j
        write ( *, '(a,2g14.6)' ) &
          'XC,YC(I+1,J) = XC,YC(I+1,J+1)=',xc(i+1,j),yc(i+1,j)
      else if ( xc(i+1,j+1) == xc(i,j+1) .and. yc(i+1,j+1) == yc(i,j+1) ) then
        nbad = nbad+1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETARE - Fatal error!'
        write ( *, '(a,2i6)' ) '  Zero length side for cell I,J:',i,j
        write ( *, '(a,2g14.6)' ) &
          'XC,YC(I+1,J+1) = XC,YC(I,J+1)=',xc(i+1,j+1),yc(i+1,j+1)
      else if ( xc(i,j+1) == xc(i,j) .and. yc(i,j+1) == yc(i,j) ) then
        nbad = nbad+1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETARE - Fatal error!'
        write ( *, '(a,2i6)' ) '  Zero length side for cell I,J:',i,j
        write ( *, '(a,2g14.6)' ) &
          'XC,YC(I,J+1) = XC,YC(I,J)=',xc(i,j+1),yc(i,j+1)
      end if
    end do
  end do

  if ( nbad > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETARE - Fatal error!'
    write ( *, '(a,i6,a)' ) '  A total of ',nbad,' zero length cell sides.'
    stop
  end if
!
!  Check for illegal zero area cells.
!
  nbad = 0
  do i = 2,l0-1
    do j = 2,m0-1
      if ( area(i,j) == 0.0D+00 ) then
        nbad = nbad+1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SETARE - Fatal error!'
        write ( *, '(a,i6,a,i6)' ) '  Zero area for cell I = ',i,' J=',j
      end if
    end do
  end do

  if ( nbad > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SETARE - Fatal error!'
    write ( *, '(a,i6,a)' ) '  A total of ',nbad,' null cells.'
    stop
  end if
!
!  Compute liquid, solid, and total areas.
!
  areal = 0.0
  do i = 1,icrys-1
    do j = 1,m0
      areal = areal+area(i,j)
    end do
  end do

  areas = 0.0
  do i = icrys,l0
    do j = 1,jcrys
      areas = areas+area(i,j)
    end do
  end do

  areat = 0.0
  do i = 1,l0
    do j = 1,m0
      areat = areat+area(i,j)
    end do
  end do

  return
end
subroutine findp ( j, xold, yold, xnew, ynew, xc, yc )

!*****************************************************************************80
!
!! FINDP finds the boundary nodes associated with a Neumann condition.
!
!  Discussion:
!
!    FINDP is given a point (XOLD,YOLD).
!
!    FINDP returns a pair of values (XNEW,YNEW), which represent ???
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  real ( kind = 8 ) a1
  real ( kind = 8 ) ai
  real ( kind = 8 ) aj
  real ( kind = 8 ) alp
  real ( kind = 8 ) b1
  real ( kind = 8 ) bet
  real ( kind = 8 ) c1
  real ( kind = 8 ) dels
  real ( kind = 8 ) denm
  real ( kind = 8 ) dxb1
  real ( kind = 8 ) dxb2
  real ( kind = 8 ) dyb1
  real ( kind = 8 ) dyb2
  integer ( kind = 4 ) j
  real ( kind = 8 ) r0
  real ( kind = 8 ) r1
  real ( kind = 8 ) r11
  real ( kind = 8 ) r2
  real ( kind = 8 ) r22
  real ( kind = 8 ) r33
  real ( kind = 8 ) slpi
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xnew
  real ( kind = 8 ) xnew1
  real ( kind = 8 ) xnew2
  real ( kind = 8 ) xold
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yn1
  real ( kind = 8 ) yn2
  real ( kind = 8 ) ynew
  real ( kind = 8 ) yold
!
!  Consider the three points P, Q and R, which have a constant
!  KSI coordinate, and an increasing ETA coordinate:
!
!    P = ( XC(2,J-1), YC(2,J-1) )
!    Q = ( XC(2,J),   YC(2,J) )
!    R = ( XC(2,J+1), YC(2,J+1) )
!
!
!  ^      R  =  [2,J+1]
!  |
!  E
!  T      Q  =  [2,J]
!  A
!  |
!  |      P  =  [2,J-1]
!  |
!  +------KSI------------->
!
!  The slope of the line from P to Q is (YQ-YP)/(XQ-XP), and
!  the slope of the line from Q to R is (YR-YQ)/(XR-XQ).
!
!  If these slopes are close enough, we can assume the three points
!  lie on a straight line.  So look at the size of
!    DENM  =  (YQ-YP)*(XR-XQ)-(YR-YQ)*(XQ-XP).
!
  dxb1 = xc(2,j)-xc(2,j-1)
  dxb2 = xc(2,j+1)-xc(2,j)
  dyb1 = yc(2,j)-yc(2,j-1)
  dyb2 = yc(2,j+1)-yc(2,j)
  denm = dyb2*dxb1-dyb1*dxb2
!
  if ( abs(denm) < 0.0001D+00 ) then

    if ( abs(dxb1) >= 0.001D+00 ) then
      slpi = -dyb1/dxb1
      aj = yc(2,j)+slpi*xc(2,j)
      ai = xold-slpi*yold
      ynew = (aj-ai*slpi)/(1.0+slpi*slpi)
      xnew = slpi*ynew+ai
    else
      ynew = yold
      xnew = xc(2,j)
    end if
!
!  Use 3 points to find a circle.
!  (X-ALP)**2+(Y-BET)**2 = R0
!
!
!  Find the cross point between circle and straight line
!  x = ai+slpi*y get a1*y**2-2*b1*y+c1=0.
!
  else

    r11 = xc(2,j-1)**2+yc(2,j-1)**2
    r22 = xc(2,j)**2+yc(2,j)**2
    r33 = xc(2,j+1)**2+yc(2,j+1)**2

    r1 = 0.5D+00*(r22-r11)
    r2 = 0.5D+00*(r33-r22)

    alp = (r1*dyb2-r2*dyb1)/denm
    bet = (r2*dxb1-r1*dxb2)/denm

    r0 = (alp-xc(2,j))**2+(bet-yc(2,j))**2

    slpi = (xold-alp)/(yold-bet)
    ai = alp-slpi*bet

    a1 = 1.0D+00+slpi*slpi
    b1 = slpi*(alp-ai)+bet
    c1 = (alp-ai)**2+bet*bet-r0

    dels = sqrt(b1*b1-a1*c1)

    yn1 = (b1+dels)/a1
    yn2 = (b1-dels)/a1

    xnew1 = slpi*yn1+ai
    xnew2 = slpi*yn2+ai

    r1 = (xnew1-xc(2,j))**2+(yn1-yc(2,j))**2
    r2 = (xnew2-xc(2,j))**2+(yn2-yc(2,j))**2

    if ( r1 < r2 ) then
      ynew = yn1
    else
      ynew = yn2
    end if

    xnew = slpi*ynew+ai

  end if

  return
end
subroutine flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
  b3jbl,cappa,cd,cmu,con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta, &
  hksi,icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho,m1,mode,nf, &
  nsolve,ntimes,r,relax,res,rueta,ruksi)

!*****************************************************************************80
!
!! FLUX computes the flux of various quantities.
!
!  It looks like if ISOL is 0, you just go through the big loop once.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: nk = 14
  integer ( kind = 4 ), parameter :: nmaxij = 64
  integer ( kind = 4 ), parameter :: ns = 10
!
  real ( kind = 8 ) acof
  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) ap0(ni,nj)
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) b3jbl(ni)
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) cmu
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) con0(ni,nj)
  real ( kind = 8 ) diff
  real ( kind = 8 ) epsil
  real ( kind = 8 ) error
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fjbi1(nj,ns)
  real ( kind = 8 ) fjbj1(ni,ns)
  real ( kind = 8 ) fjbl1(nj,ns)
  real ( kind = 8 ) fjbm1(ni,ns)
  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) flow
  real ( kind = 8 ) fmax(ns)
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) isol
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) izone
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) l1
  logical lblk(ns)
  logical lconv
  logical lortho
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nsolve(ns)
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) ntimes(ns)
  real ( kind = 8 ) pt(nmaxij)
  real ( kind = 8 ) qt(nmaxij)
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) tem1
  real ( kind = 8 ) tem2
  real ( kind = 8 ) xje
  real ( kind = 8 ) xjn
  real ( kind = 8 ) xjw
  real ( kind = 8 ) xjs
  real ( kind = 8 ) xme
  real ( kind = 8 ) xmn
  real ( kind = 8 ) xmw
  real ( kind = 8 ) xms
!
!  Save copies of CON and AP.
!
  do i = 1,l1
    do j = 1,m1
      con0(i,j) = con(i,j)
      ap0(i,j) = ap(i,j)
    end do
  end do
!
!  I think ITER = 0 occurs only for Temperature in the Solid zone...
!
  if ( iter == 0 ) then
    if ( isol /= 0 ) then

      call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
        heta,hksi,icrys,izone,jcrys,l1,m1,nf)

      call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)

    end if
    return
  end if

  do nt = 1,ntimes(nf)
!
!  Find the flux on the first step only.
!
    if ( nt == 1 ) then

      do j = 1,m1

        diff = gam(2,j)/(0.5D+00*hksi(2,j))
        flow = -ruksi(2,j)
        call diflow(acof,diff,flow)
        qt(2) = ruksi(2,j)*f(2,j,nf)+acof*(f(1,j,nf)-f(2,j,nf))

        do i = 3,l1-1
          jj = j
          if ( j == 1) jj = 2
          if ( j == m1) jj = m1-1
          diff = gam(i,jj)*gam(i-1,jj)/(0.5D+00*hksi(i,jj)*gam(i-1,jj) &
            +0.5D+00*hksi(i-1,jj)*gam(i,jj))
          flow = -ruksi(i,j)
          call diflow(acof,diff,flow)
          qt(i) = ruksi(i,j)*f(i,j,nf)+acof*(f(i-1,j,nf)-f(i,j,nf))
        end do

        diff = gam(l1-1,j)/(0.5D+00*hksi(l1-1,j))
        flow = ruksi(l1,j)
        call diflow(acof,diff,flow)
        qt(l1) = ruksi(l1,j)*f(l1-1,j,nf)+acof*(f(l1-1,j,nf)-f(l1,j,nf))

        do i = 2,l1
          fjksi(i,j) = qt(i)
        end do

      end do
!
!  Along J.
!
      do i = 1,l1

        diff = gam(i,2)/(0.5D+00*heta(i,2))
        flow = -rueta(i,2)
        call diflow(acof,diff,flow)
        qt(2) = rueta(i,2)*f(i,2,nf)+acof*(f(i,1,nf)-f(i,2,nf))

        do j = 3,m1-1

          if ( i == 1 ) then
            ii = 2
          else if ( i == l1 ) then
            ii = l1-1
          else
            ii = i
          end if

          diff = gam(ii,j)*gam(ii,j-1)/(0.5D+00*heta(ii,j)*gam(ii,j-1) &
            +0.5D+00*heta(ii,j-1)*gam(ii,j))
          flow = -rueta(i,j)
          call diflow(acof,diff,flow)
          qt(j) = rueta(i,j)*f(i,j,nf)+acof*(f(i,j-1,nf)-f(i,j,nf))

        end do

        diff = gam(i,m1-1)/(0.5D+00*heta(i,m1-1))
        flow = rueta(i,m1)
        call diflow(acof,diff,flow)
        qt(m1) = rueta(i,m1)*f(i,m1-1,nf)+acof*(f(i,m1-1,nf)-f(i,m1,nf))

        do j = 2,m1
          fjeta(i,j) = qt(j)
        end do

      end do
!
!  Finish JKSI and JETA calculation.
!  Construct boundary flux for J type boundary condition.
!  Boundary condition con(1,j),*,*,* go into jksi and jeta
!  ak2/ak1,* is project on xi - direction from eta - direction
!  variable, in program let ak2 = 0.0, because sym. orth.
!
      do j = 2,m1-1

        if ( gam(1,j) == 0.0D+00 ) then
          fjksi(2,j) = (con(1,j)*heta(1,j)+0.5D+00*ak2(2,j) &
            *(fjeta(1,j)+fjeta(1,j+1)))/ak1(2,j)
        end if

        if ( gam(l1,j) == 0.0D+00 ) then
          fjksi(l1,j) = (con(l1,j)*heta(l1,j)+0.5D+00*ak2(l1,j) &
            *(fjeta(l1,j)+fjeta(l1,j+1)))/ak1(l1,j)
        end if

      end do

      do i = 2,l1-1

        if ( gam(i,1) == 0.0D+00 ) then
          fjeta(i,2) = (con(i,1)*hksi(i,1)+0.5D+00*ae2(i,2) &
            *(fjksi(i,1)+fjksi(i+1,1)))/ae1(i,2)
        end if

        if ( gam(i,m1) == 0.0D+00 ) then
          fjeta(i,m1) = (con(i,m1)*hksi(i,m1)+0.5D+00*ae2(i,m1) &
            *(fjksi(i,m1)+fjksi(i+1,m1)))/ae1(i,m1)
        end if

      end do
!
!  (2-31) and (2-42) find out vector J cdot vector n or J_n in bound.
!
      do i = 2,l1-1

        t1 = 0.5D+00*(fjksi(i,1)+fjksi(i+1,1))
        t2 = 0.5D+00*(fjksi(i,m1)+fjksi(i+1,m1))
        fjbj1(i,nf) = (fjeta(i,2)*ae1(i,2)-t1*ae2(i,2))/hksi(i,1)
        fjbm1(i,nf) = (fjeta(i,m1)*ae1(i,m1)-t2*ae2(i,m1))/hksi(i,m1)

        if ( mode == 1 ) then
          fjbj1(i,nf) = fjbj1(i,nf)/r(i,2)
          fjbm1(i,nf) = fjbm1(i,nf)/r(i,m1)
        end if

      end do

      do j = 2,m1-1

        t1 = 0.5D+00*(fjeta(1,j)+fjeta(1,j+1))
        t2 = 0.5D+00*(fjeta(l1,j)+fjeta(l1,j+1))
        fjbi1(j,nf) = (fjksi(2,j)*ak1(2,j)-t1*ak2(2,j))/heta(1,j)
        fjbl1(j,nf) = (fjksi(l1,j)*ak1(l1,j)-t2*ak2(l1,j))/heta(l1,j)

        if ( mode == 1 ) then
          fjbi1(j,nf) = fjbi1(j,nf)/r(2,j)
          fjbl1(j,nf) = fjbl1(j,nf)/r(l1,j)
        end if

      end do
!
!  This code is only for the temperature variable.
!  We compute
!    B1JBL, the heat flux between rows (?) ICRYS-1 and ICRYS,
!           for columns 2 through JCRYS and column M1.
!    B2JBL, the heat flux between rows (?) ICRYS and ICRYS+1,
!           for columns 2 through JCRYS and column M1.
!    B3JBL, the heat flux between rows (?) and (?),
!           for columns JCRYS+1 to M1-1.
!
      if ( nf == 5 ) then

        do j = 2,jcrys-1

          t1 = gam(icrys-1,j)*(f(icrys-1,j,5)-f(icrys,j,5))/hksi(icrys-1,j)

          b1jbl(j) = t1*ak1(icrys,j)/heta(icrys,j)

          if ( mode == 1 ) then
            b1jbl(j) = b1jbl(j)/r(icrys,j)
          end if

        end do

        b1jbl(jcrys) = b1jbl(jcrys+1)
        b1jbl(m1) = b1jbl(m1-1)

        do j = 2,jcrys-1

          t1 = gam(icrys+1,j)*(f(icrys,j,5)-f(icrys+1,j,5))/hksi(icrys+1,j)

          b2jbl(j) = t1*ak1(icrys,j)/heta(icrys,j)

          if ( mode == 1 ) then
            b2jbl(j) = b2jbl(j)/r(icrys,j)
          end if

        end do

        b2jbl(jcrys) = b2jbl(jcrys+1)
        b2jbl(m1) = b2jbl(m1-1)

        do j = jcrys+1,m1-1

          t1 = gam(icrys-1,j)*(f(icrys-1,j,5)-f(icrys,j,5))/hksi(icrys-1,j)

          t2 = gam(icrys-1,j)*(f(icrys,j+1,5)-f(icrys,j,5))/heta(icrys,j)

          b3jbl(j) = (t1*ak1(icrys,j)-t2*ak2(icrys,j))/heta(icrys,j)

          if ( mode == 1 ) then
            b3jbl(j) = b3jbl(j)/r(icrys,j)
          end if

        end do

        b3jbl(jcrys) = b3jbl(jcrys+1)
        b3jbl(m1) = b3jbl(m1-1)

      end if

    else
!
!  This code is done if this is NOT the first iteration.
!
!  Correct J from known PHI and J.
!  Correct JKSI.
!
      do j = 2,m1-1

        if ( gam(1,j) /= 0.0D+00 ) then
          fjksi(2,j) = fjksi(2,j)+aim(2,j)/ak1(2,j)* &
            (f(1,j,nf)-f(2,j,nf))+ruksi(2,j)*f(2,j,nf)
        end if

        do i = 3,l1

          if ( gam(l1,j) /= 0.0D+00 .or. i.ne.l1 ) then
            fjksi(i,j) = fjksi(i,j)+aip(i-1,j)/ak1(i,j)*(f(i-1,j,nf) &
              -f(i,j,nf))+ruksi(i,j)*f(i-1,j,nf)
          end if

        end do

      end do
!
!  Correct FJETA.
!
      do i = 2,l1-1

        if ( gam(i,1) /= 0.0D+00 ) then
          fjeta(i,2) = fjeta(i,2)+ajm(i,2)/ae1(i,2)* &
            (f(i,1,nf)-f(i,2,nf))+rueta(i,2)*f(i,2,nf)
        end if

        do j = 3,m1

          if ( gam(i,m1) /= 0.0D+00 .or. j.ne.m1 ) then
            fjeta(i,j) = fjeta(i,j)+ajp(i,j-1)/ae1(i,j)* &
              (f(i,j-1,nf)-f(i,j,nf))+rueta(i,j)*f(i,j-1,nf)
          end if

        end do
      end do

    end if
!
!  End of "Is this the first iteration or not" block.
!
!  Restore old values of CON and AP.
!
    do i = 1,l1
      do j = 1,m1
        con(i,j) = con0(i,j)
        ap(i,j) = ap0(i,j)
      end do
    end do

    if ( .not.lortho ) then

      do j = 2,m1-1

        do i = 1,l1
          pt(i) = 0.5D+00*(fjeta(i,j)+fjeta(i,j+1))
          qt(i) = 0.5D+00*(rueta(i,j)+rueta(i,j+1))
        end do

        do i = 2,l1-1
          t1 = hksi(i,j)/(hksi(i+1,j)+hksi(i,j))
          t2 = hksi(i+1,j)/(hksi(i+1,j)+hksi(i,j))
          t3 = hksi(i,j)/(hksi(i-1,j)+hksi(i,j))
          t4 = hksi(i-1,j)/(hksi(i-1,j)+hksi(i,j))
          xje = t2*pt(i)+t1*pt(i+1)
          xme = t2*qt(i)+t1*qt(i+1)
          xjw = t4*pt(i)+t3*pt(i-1)
          xmw = t4*qt(i)+t3*qt(i-1)
          con(i,j) = con(i,j)+(xje-f(i,j,nf)*xme)*ak2(i+1,j) &
            -(xjw-f(i,j,nf)*xmw)*ak2(i,j)
        end do
      end do

      do i = 2,l1-1

        do j = 1,m1
          pt(j) = 0.5D+00*(fjksi(i,j)+fjksi(i+1,j))
          qt(j) = 0.5D+00*(ruksi(i,j)+ruksi(i+1,j))
        end do

        do j = 2,m1-1
          tem1 = heta(i,j+1)+heta(i,j)
          tem2 = heta(i,j-1)+heta(i,j)
          t1 = heta(i,j)/tem1
          t2 = heta(i,j+1)/tem1
          t3 = heta(i,j)/tem2
          t4 = heta(i,j-1)/tem2
          xjn = t2*pt(j)+t1*pt(j+1)
          xmn = t2*qt(j)+t1*qt(j+1)
          xjs = t4*pt(j)+t3*pt(j-1)
          xms = t4*qt(j)+t3*qt(j-1)
          con(i,j) = con(i,j)+(xjn-f(i,j,nf)*xmn)*ae2(i,j+1) &
            -(xjs-f(i,j,nf)*xms)*ae2(i,j)
        end do
      end do

    end if
!
!  Calculate JHAT.
!  Step 5 use equation (2-97) and (2-102), (2-103), (2-104)
!  calculate boundary data from J and other purpose for B_SP.
!  from here to end. pc  =  JKSI HAT or JETA HAT.
!
!  jksi hat
!
    do j = 2,m1-1

      f(2,j,4) = 0.0D+00

      if ( gam(1,j) == 0.0D+00 ) then

        t1 = fjksi(2,j)-ruksi(2,j)*f(2,j,nf)
        diff = gam(2,j)/(0.5D+00*hksi(2,j))
        flow = -ruksi(2,j)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(1,j,nf) = f(2,j,nf)+t1/acof
        else
          f(1,j,nf) = f(2,j,nf)
        end if

        f(2,j,4) = fjksi(2,j)

      end if

      f(l1,j,4) = 0.0D+00

      if ( gam(l1,j) == 0.0D+00 ) then

        t1 = fjksi(l1,j)-ruksi(l1,j)*f(l1-1,j,nf)
        diff = gam(l1-1,j)/(0.5D+00*hksi(l1-1,j))
        flow = ruksi(l1,j)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(l1,j,nf) = f(l1-1,j,nf)-t1/acof
        else
          f(l1,j,nf) = f(l1-1,j,nf)
        end if

        f(l1,j,4) = fjksi(l1,j)

      end if

      do i = 3,l1-1
        f(i,j,4) = 0.0D+00
      end do

    end do

    do j = 2,m1-1
      do i = 2,l1
        fjksi(i,j) = f(i,j,4)
      end do
    end do
!
!  jeta hat
!
    do i = 2,l1-1

      f(i,2,4) = 0.0D+00

      if ( gam(i,1) == 0.0D+00 ) then

        t1 = fjeta(i,2)-rueta(i,2)*f(i,2,nf)
        diff = gam(i,2)/(0.5D+00*heta(i,2))
        flow = -rueta(i,2)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(i,1,nf) = f(i,2,nf)+t1/acof
        else
          f(i,1,nf) = f(i,2,nf)
        end if

        f(i,2,4) = fjeta(i,2)

      end if

      f(i,m1,4) = 0.0D+00

      if ( gam(i,m1) == 0.0D+00 ) then

        t1 = fjeta(i,m1)-rueta(i,m1)*f(i,m1-1,nf)
        diff = gam(i,m1-1)/(0.5D+00*heta(i,m1-1))
        flow = rueta(i,m1)
        call diflow(acof,diff,flow)

        if ( acof /= 0.0D+00 ) then
          f(i,m1,nf) = f(i,m1-1,nf)-t1/acof
        else
          f(i,m1,nf) = f(i,m1-1,nf)
        end if

        f(i,m1,4) = fjeta(i,m1)

      end if

      do j = 3,m1-1
        f(i,j,4) = 0.0D+00
      end do

    end do

    do i = 2,l1-1
      do j = 2,m1
        fjeta(i,j) = f(i,j,4)
      end do
    end do
!
!  Calculate and solve phy equation
!  (2-106) and (2-107) b  =  b_s + b_no + b_sp, t1=b_sp
!
    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = con(i,j)+fjksi(i,j)*ak1(i,j)-fjksi(i+1,j)*ak1(i+1,j) &
          +fjeta(i,j)*ae1(i,j)-fjeta(i,j+1)*ae1(i,j+1)
      end do
    end do

    if ( isol == 0) return

    do i = 2,l1-1
      do j = 2,m1-1
        ap(i,j) = ap(i,j)/relax(nf)
        con(i,j) = con(i,j)+(1.0D+00-relax(nf))*ap(i,j)*f(i,j,nf)
      end do
    end do

    call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
      heta,hksi,icrys,izone,jcrys,l1,m1,nf)

    call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)
!
!  BEGIN NEW.
!
!  Compute the largest solution component.
!
    fmax(nf) = 0.0
    do i = 1,l1
      do j = 1,m1
        fmax(nf) = max(fmax(nf),abs(f(i,j,nf)))
      end do
    end do
!
!  Compare current and previous iterates.
!
    if ( nt > 1 ) then
      res(nf) = 0.0D+00
      do i = 1,l1
        do j = 1,m1

          error = abs(f(i,j,nf)-fn(i,j,nf))

          if ( fmax(nf) > 0.0D+00 ) then
            error = error/fmax(nf)
          end if

          res(nf) = max(res(nf),error)

        end do
      end do
    end if
!
!  Save the new solution in FN.
!
    do i = 1,l1
      do j = 1,m1
        fn(i,j,nf) = f(i,j,nf)
      end do
    end do

    if ( res(nf) < epsil .and. nt > 1 ) then
!
!         write ( *, '(a,i6)' ) 'FLUX - Convergence for variable ',nf
!         write ( *, '(a,i6)' ) '  on step ',nt
!         write ( *, '(a,g14.6)' ) '  RES(NF) = ',res(nf)
!         write ( *, '(a,g14.6)' ) '  FMAX(NF)=',fmax(nf)
!
      return
    end if

  end do
!
!  End of big iteration
!
!  Compute the largest solution component.
!
!     fmax(nf) = 0.0D+00
!     do i = 1,l1
!       do j = 1,m1
!         fmax(nf) = max(fmax(nf),abs(f(i,j,nf)))
!       end do
!     end do
!
!  Compute the largest relative change in the solution.
!
!     res(nf) = 0.0D+00
!     do i = 1,l1
!       do j = 1,m1
!
!         error = abs(f(i,j,nf)-fn(i,j,nf))
!
!         if ( fmax(nf) > 0.0D+00 ) then
!           error = error/fmax(nf)
!         end if
!
!         res(nf) = max(res(nf),error)
!
!       end do
!     end do
!
!  Decide whether to declare convergence or not.
!
  if ( res(nf) > epsil ) then
    lconv = .false.
  end if
!
!  Save the new solution in FN.
!
!     do i = 1,l1
!       do j = 1,m1
!         fn(i,j,nf) = f(i,j,nf)
!       end do
!     end do

!     nt = ntimes(nf)
!     write ( *, '(a,i6)' ) 'FLUX - No convergence for variable ',nf
!     write ( *, '(a,i6)' ) '  on step ',nt
!     write ( *, '(a,g14.6)' ) '  RES(NF) = ',res(nf)
!     write ( *, '(a,g14.6)' ) '  FMAX(NF)=',fmax(nf)

  return
end
subroutine gamsor(ap,birad,ce1,ce2,cmu,con,ewall,f,fcsl,fksl, &
  fma,fo,frsl,gam,gamt,grash,hamag,heta,hksi,icrys,inturb, &
  izone,jcrys,l1,lsolve,m1,mode,nf,pr,r,rdtm,re,recb,rect, &
  rho,rhocon,rpr,sige,sigk,sigt,stel,stes,tal,tas,tf,tw,x, &
  xc,y,yc)
!
!*****************************************************************************80
!
!! GAMSOR sets the coefficients for the transport problems.
!
!  Discussion:
!
!    These coefficients are associated with U, V, T, W, TK, TE and E.
!
!    Note that GAM(I,1) = 0 means the gradient is zero at wall J=1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: nmaxij = 64
  integer ( kind = 4 ), parameter :: ns = 10

  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) ar
  real ( kind = 8 ) birad
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cm4
  real ( kind = 8 ) cmu
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) dcdr
  real ( kind = 8 ) depdr
  real ( kind = 8 ) depdx
  real ( kind = 8 ) dudx(ni,nj)
  real ( kind = 8 ) dudy(ni,nj)
  real ( kind = 8 ) dvdx(ni,nj)
  real ( kind = 8 ) dvdy(ni,nj)
  real ( kind = 8 ) dwdx(ni,nj)
  real ( kind = 8 ) dwdy(ni,nj)
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) frsl
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) gdudx(ni,nj)
  real ( kind = 8 ) gdudy(ni,nj)
  real ( kind = 8 ) gdvdx(ni,nj)
  real ( kind = 8 ) gdvdy(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) inturb
  integer ( kind = 4 ) izone
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l1
  logical lsolve(ns)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) modelm
  integer ( kind = 4 ) nf
  real ( kind = 8 ) p11
  real ( kind = 8 ) p12
  real ( kind = 8 ) p21
  real ( kind = 8 ) p22
  real ( kind = 8 ) p31
  real ( kind = 8 ) p32
  real ( kind = 8 ) p41
  real ( kind = 8 ) p42
  real ( kind = 8 ) pbig
  real ( kind = 8 ) pr
  real ( kind = 8 ) prod(ni,nj)
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) rpr
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigt
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) tal
  real ( kind = 8 ) tas
  real ( kind = 8 ) tf
  real ( kind = 8 ) tw
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xplus(nmaxij)
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yplus(nmaxij)
!
  do i = 1,l1
    do j = 1,m1
      ap(i,j) = -rdtm*rho(i,j)
      con(i,j) = 0.0D+00
      gam(i,j) = rhocon/re+gamt(i,j)
    end do
  end do
!
!  Horizontal velocity.
!
  if ( nf == 1 ) then

    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = fo(i,j,1)*rho(i,j)*rdtm+rho(i,j)*f(i,j,5)*grash/re**2
      end do
    end do

    if ( inturb == 1 ) then

      do i = 1,l1
        do j = 1,m1
          gdudx(i,j) = 0.0D+00
          gdudy(i,j) = 0.0D+00
          gdvdx(i,j) = 0.0D+00
          gdvdy(i,j) = 0.0D+00
        end do
      end do

      do i = 2,l1-1
        do j = 2,m1-1
          call gradnt(heta,hksi,i,j,f(1,1,1),dudx(i,j),dudy(i,j),xc,yc)
          call gradnt(heta,hksi,i,j,f(1,1,2),dvdx(i,j),dvdy(i,j),xc,yc)
          gdudx(i,j) = gamt(i,j)*dudx(i,j)
          gdudy(i,j) = gamt(i,j)*dudy(i,j)
          gdvdx(i,j) = gamt(i,j)*dvdx(i,j)
          gdvdy(i,j) = gamt(i,j)*dvdy(i,j)
        end do
      end do

      do i = 2,l1-1
        do j = 2,m1-1
          call gradnt(heta,hksi,i,j,gdudx,p11,p12,xc,yc)
          call gradnt(heta,hksi,i,j,gdudy,p21,p22,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdx,p31,p32,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdy,p41,p42,xc,yc)
          con(i,j) = con(i,j)+p11+p32
        end do
      end do

      cm4 = cmu**0.25

      if ( mode == 0 ) then
        do i = 2,l1-1

          yplus(i) = re*rho(i,2)*sqrt(f(i,2,7))*cm4*0.5*heta(i,2)/rhocon

          if ( yplus(i) <= 11.63D+00 ) then
            gam(i,1) = rhocon/re
          else
            gam(i,1) = yplus(i)/(2.5D+00*log(ewall*yplus(i)))
          end if

        end do
      end if

      do i = 2,l1-1

        yplus(i) = re*rho(i,m1-1)*sqrt(f(i,m1-1,7)) &
          *cm4*0.5D+00*heta(i,m1-1)/rhocon

        if ( yplus(i) <= 11.63D+00 ) then
          gam(i,m1) = rhocon/re
        else
          gam(i,m1) = yplus(i)/(2.5D+00*log(ewall*yplus(i)))
        end if

      end do

      do j = 2,m1-1

        xplus(j) = re*rho(2,j)*sqrt(f(2,j,7))*cm4*0.5D+00*hksi(2,j)/rhocon

        if ( xplus(j) <= 11.63D+00 ) then
          gam(1,j) = rhocon/re
        else
          gam(1,j) = xplus(j)/(2.5D+00*log(ewall*xplus(j)))
        end if

      end do

      do j = 2,m1-1

        xplus(j) = re*rho(l1-1,j)*sqrt(f(l1-1,j,7)) &
          *cm4*0.5D+00*hksi(l1-1,j)/rhocon

        if ( xplus(j) <= 11.63D+00 ) then
          gam(l1,j) = rhocon/re
        else
          gam(l1,j) = xplus(j)/(2.5D+00*log(ewall*xplus(j)))
        end if

      end do

    end if

    do i = 2,l1-1
      gam(i,1) = 0.0D+00
    end do

    do j = 2,m1-1
      f(1,j,1) = 0.0D+00
      f(l1,j,1) = 0.0D+00
    end do

    f(l1,1,1) = 0.0D+00
    do j = jcrys+1,m1-1
      f(l1,j,1) = fo(l1,j,2)*(xc(l1,j+1)-xc(l1,j))/(yc(l1,j+1)-yc(l1,j))
    end do
    f(l1,jcrys+1,1) = 0.0D+00
!
!  Vertical velocity.
!
  else if ( nf == 2 ) then

    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = fo(i,j,2)*rho(i,j)*rdtm
        ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re
        if ( mode == 1 ) then
          con(i,j) = con(i,j)+rho(i,j)*f(i,j,6)**2/r(i,j)**3
          ap(i,j) = ap(i,j)-rho(i,j)/r(i,j)**2
        end if
      end do
    end do

    if ( inturb == 1 ) then

      do i = 1,l1
        do j = 1,m1
          gdudx(i,j) = 0.0D+00
          gdudy(i,j) = 0.0D+00
          gdvdx(i,j) = 0.0D+00
          gdvdy(i,j) = 0.0D+00
        end do
      end do

      do i = 2,l1-1
        do j = 2,m1-1
          call gradnt(heta,hksi,i,j,f(1,1,1),dudx(i,j),dudy(i,j),xc,yc)
          call gradnt(heta,hksi,i,j,f(1,1,2),dvdx(i,j),dvdy(i,j),xc,yc)
          gdudx(i,j) = gamt(i,j)*dudx(i,j)
          gdudy(i,j) = gamt(i,j)*dudy(i,j)
          gdvdx(i,j) = gamt(i,j)*dvdx(i,j)
          gdvdy(i,j) = gamt(i,j)*dvdy(i,j)
        end do
      end do

      do i = 2,l1-1
        do j = 2,m1-1
          call gradnt(heta,hksi,i,j,gdudx,p11,p12,xc,yc)
          call gradnt(heta,hksi,i,j,gdudy,p21,p22,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdx,p31,p32,xc,yc)
          call gradnt(heta,hksi,i,j,gdvdy,p41,p42,xc,yc)
          con(i,j) = con(i,j)+p21+p42
        end do
      end do

      cm4 = cmu**0.25

      if ( mode == 0 ) then
        do i = 2,l1-1

          yplus(i) = re*rho(i,2)*sqrt(f(i,2,7))*cm4*0.5D+00*heta(i,2)/rhocon
          gam(i,1) = rhocon/re
          if ( yplus(i) > 11.63D+00) then
            gam(i,1) = yplus(i)/(2.5D+00*log(ewall*yplus(i)))
          end if

        end do
      end if

      do i = 2,l1-1

        yplus(i) = re*rho(i,m1-1)*sqrt(f(i,m1-1,7)) &
          *cm4*0.5*heta(i,m1-1)/rhocon
        gam(i,m1) = rhocon/re
        if ( yplus(i) > 11.63D+00 ) then
          gam(i,m1) = yplus(i)/(2.5D+00*log(ewall*yplus(i)))
        end if

      end do

      do j = 2,m1-1
        xplus(j) = re*rho(2,j)*sqrt(f(2,j,7))*cm4*0.5D+00*hksi(2,j)/rhocon
        gam(1,j) = rhocon/re
        if ( xplus(j) > 11.63D+00 ) then
          gam(1,j) = xplus(j)/(2.5D+00*log(ewall*xplus(j)))
        end if
      end do

      do j = 2,m1-1
        xplus(j) = re*rho(l1-1,j)*sqrt(f(l1-1,j,7)) &
          *cm4*0.5D+00*hksi(l1-1,j)/rhocon
        gam(l1,j) = rhocon/re
        if ( xplus(j) > 11.63D+00 ) then
          gam(l1,j) = xplus(j)/(2.5D+00*log(ewall*xplus(j)))
        end if
      end do

    end if

    do i = 2,l1-1
      f(i,1,2) = 0.0D+00
    end do

    do j = jcrys+1,m1-1
      dcdr = (f(l1,j+1,5)-f(l1,j,5))/(y(l1,j+1)-y(l1,j))
      f(l1,j,2) = f(l1-1,j,2)+(xc(l1,j)-xc(l1-1,j))*dcdr*(fma/(re*pr))
    end do

    f(l1,jcrys+1,2) = 0.0D+00
    do j = 2,jcrys
      f(1,j,2) = 0.0D+00
      f(l1,j,2) = 0.0D+00
    end do

    f(l1,m1-1,2) = 0.0D+00
    f(l1,m1,2) = 0.0D+00
!
!  Temperature computations in the solid zone.
!
  else if ( nf == 5 .and. izone == 1 ) then

    do i = 1,l1
      do j = 1,m1
        gam(i,j) = rpr*fksl/fcsl/frsl
      end do
    end do

    do i = icrys+1,l1
      gam(i,1) = 0.0D+00
    end do

    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = fo(i,j,5)*rho(i,j)*rdtm
      end do
    end do

    do j = 2,jcrys
      f(l1,j,5) = -stes/stel
    end do

    do i = icrys+1,l1
      f(i,m1,5) = f(i,m1-1,5)-birad*((f(i,m1,5)*(tw-tf)+tf)**4-tas**4) &
        *(y(i,m1)-y(i,m1-1))/100000000.0D+00
    end do
!
!  Temperature calculations in the liquid zone.
!
  else if ( nf == 5 .and. izone == 2 ) then

    do i = 1,l1
      do j = 1,m1
        gam(i,j) = rpr+gamt(i,j)/sigt
      end do
    end do

    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = fo(i,j,5)*rho(i,j)*rdtm
      end do
    end do

    if ( inturb == 1 ) then

      cm4 = cmu**0.25
      pbig = 9.24D+00*((pr/sigt)**0.75-1.0D+00) &
        *(1.0D+00+0.28D+00*exp(-0.007D+00*pr/sigt))

      if ( mode == 0) then
        do i = 2,l1-1
          yplus(i) = re*rho(i,2)*sqrt(f(i,2,7))*cm4*0.5D+00*heta(i,2)/rhocon
          gam(i,1) = rpr
          if ( yplus(i) > 11.63D+00 ) then
            gam(i,1) = max(rhocon/re, &
              yplus(i)/(sigt*(2.5D+00*log(ewall*yplus(i))+pbig)))
          end if
        end do
      end if

      do i = 2,l1-1
        yplus(i) = re*rho(i,m1-1)*sqrt(f(i,m1-1,7)) &
          *cm4*0.5D+00*heta(i,m1-1)/rhocon
        gam(i,m1) = rpr
        if ( yplus(i) > 11.63D+00 ) then
          gam(i,m1) = max(rhocon/re, &
            yplus(i)/(sigt*(2.5D+00*log(ewall*yplus(i))+pbig)))
        end if
      end do

      do j = 2,m1-1
        xplus(j) = re*rho(2,j)*sqrt(f(2,j,7))*cm4*0.5D+00*hksi(2,j)/rhocon
        gam(1,j) = rpr
        if ( xplus(j) > 11.63D+00 ) then
          gam(1,j) = max(rhocon/re, &
            xplus(j)/(sigt*(2.5*log(ewall*xplus(j))+pbig)))
        end if
      end do

      do j = 2,m1-1
        xplus(j) = re*rho(l1-1,j)*sqrt(f(l1-1,j,7)) &
          *cm4*0.5*hksi(l1-1,j)/rhocon
        gam(l1,j) = rpr
        if ( xplus(j) > 11.63D+00 ) then
          gam(l1,j) = max(rhocon/re, &
            xplus(j)/(sigt*(2.5D+00*log(ewall*xplus(j))+pbig)))
        end if
      end do

    end if

    do j = 2,jcrys
      f(l1,j,5) = 0.0D+00
    end do

    do j = jcrys+1,m1
      f(l1,j,5) = f(l1-1,j,5)-0.1D+00*birad*((f(l1,j,5)*(tw-tf)+tf)**4-tal**4) &
        *(x(l1,j)-x(l1-1,j))/100000000.0D+00
    end do

    do j = 1,m1
      f(1,j,5) = 0.5D+00 + 0.5D+00*yc(2,j)
    end do

    do i = 1,l1
      f(i,m1,5) = 1.0D+00
    end do

    do i = 1,l1
      gam(i,1) = 0.0D+00
    end do
!
!  Angular momentum.
!
  else if ( nf == 6 ) then
!
!  modelm = 1. based on rectangular de/dx of Sabhapathy and Salcudean
!  modelm = 2. based on curvilinear de/dx of Sabhapathy and Salcudean.
!  modelm = 3. same as 2 but different numerical treatment for source
!    terms based on strong magnetic. 1,2 is good for weak Ha
!    3 is good for moderate.
!  modelm = 4. based on strong magnetic. rotation is eliminated.
!  modelm = 5. is a test.
!
    if ( hamag == 0.0D+00 ) then
      modelm = 0.0D+00
    else if ( hamag > 0.0D+00 .and. hamag <= 20.0D+00 ) then
      modelm = 2
    else if ( hamag > 20.0D+00 .and. hamag <= 40.0D+00 ) then
      modelm = 5
    else if ( hamag > 40.0D+00 ) then
      modelm = 4
    end if

    do i = 2,l1-1
      do j = 2,m1-1

        con(i,j) = fo(i,j,6)*rho(i,j)*rdtm

        if ( modelm == 1 ) then

          con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)* &
            (f(i+1,j,9)-f(i,j,9))/(x(i+1,j)-x(i,j))

        else if ( modelm == 2 ) then

          call gradnt(heta,hksi,i,j,f(1,1,9),depdx,depdr,xc,yc)
          if ( j >= 3 ) then
            con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)*depdx
          end if

        else if ( modelm == 3 ) then

          call gradnt(heta,hksi,i,j,f(1,1,9),depdx,depdr,xc,yc)
          if ( j >= 3 ) then
            con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)*depdx &
              +rho(i,j)*hamag**2/re*f(i,j,6)
            ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re
          end if

        else if ( modelm == 4 ) then

          ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re

        else if ( modelm == 5 ) then

          call gradnt(heta,hksi,i,j,f(1,1,9),depdx,depdr,xc,yc)

          if ( j >= 20 ) then
            con(i,j) = con(i,j)-rho(i,j)*hamag**2/re/r(i,j)*depdx
          end if

          if ( j <= 20 ) then
            ap(i,j) = ap(i,j)-rho(i,j)*hamag**2/re
          end if

        end if
!
!  -2/r*dOmega/dr /re
!
        ar = 2.0D+00/r(i,j)/(y(i,j)-y(i,j-1))/re
        con(i,j) = con(i,j)+rho(i,j)*ar*f(i,j-1,6)
        ap(i,j) = ap(i,j)-rho(i,j)*ar

      end do
    end do

    do i = 1,l1
      f(i,m1,6) = recb
      f(i,1,6) = 0.0D+00
    end do

    do j = jcrys+1,m1
      gam(l1,j) = 0.0D+00
    end do

    do j = 1,m1

      f(1,j,6) = recb*r(2,j)**2

      if ( j > jcrys ) then
        f(l1,j,6) = f(l1-1,j,6)
      else
        f(l1,j,6) = rect*r(l1,j)**2
      end if

    end do
!
!  Turbulent kinetic energy.
!
  else if ( nf == 7 ) then

    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = fo(i,j,7)*rdtm
        gam(i,j) = rhocon/re+gamt(i,j)/sigk
      end do
    end do

    do i = 2,l1-1
      do j = 2,m1-1

        call gradnt(heta,hksi,i,j,f(1,1,1),dudx(i,j),dudy(i,j),xc,yc)
        call gradnt(heta,hksi,i,j,f(1,1,2),dvdx(i,j),dvdy(i,j),xc,yc)
        prod(i,j) = gamt(i,j)*(2.0D+00*(dudx(i,j)**2+dvdy(i,j)**2) &
          +(dudy(i,j)+dvdx(i,j))**2)

        if ( mode == 1 ) then
          prod(i,j) = prod(i,j)+gamt(i,j)*2.0D+00*(f(i,j,2)/r(i,j))**2
        end if

        if ( lsolve(6) ) then
          call gradnt(heta,hksi,i,j,f(1,1,6),dwdx(i,j),dwdy(i,j),xc,yc)
          prod(i,j) = prod(i,j)+gamt(i,j)*(dwdx(i,j)**2 &
            +(dwdy(i,j)-2.0D+00*f(i,j,6)/r(i,j))**2)/r(i,j)**2
        end if

        con(i,j) = con(i,j)+prod(i,j)
        ap(i,j) = ap(i,j)-cmu*rho(i,j)**2*abs(f(i,j,7))/(gamt(i,j)+1.0D-05)

      end do
    end do
!
!  Turbulent dissipation.
!
  else if ( nf == 8 ) then

    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = fo(i,j,8)*rdtm
        gam(i,j) = rhocon/re+gamt(i,j)/sige
      end do
    end do

    do i = 2,l1-1
      do j = 2,m1-1
        con(i,j) = con(i,j)+ce1*cmu*rho(i,j)*abs(f(i,j,7))*prod(i,j) &
          /((gamt(i,j)+1.0D-05))
        ap(i,j) = ap(i,j)-ce2*cmu*rho(i,j)**2*abs(f(i,j,7))/(gamt(i,j)+1.0D-05)
      end do
    end do

  end if

  return
end
subroutine gradnt(heta,hksi,i,j,phi,dphidx,dphidy,xc,yc)

!*****************************************************************************80
!
!! GRADNT estimats the gradient of PHI at a primary node index (I,J).
!
!  Discussion:
!
!    The routine is given:
!
!    PHI, the value of a state quantity at the primary nodes;
!    HETA, HKSI, the two dimensions of the cell, as measured
!      in the X, Y coordinate system;
!    I, J, the indices of the primary node at which a gradient
!      is desired;
!    XC, YC, the X, Y coordinates of the corner nodes that define the
!      shape of the cells, and the location of the midside nodes, and
!      primary nodes.
!
!    GRADNT must estimate
!
!    DPHIDX, DPHIDY, an estimate for the gradient (d PHI/d X, d PHI/d Y)
!      of the physical quantity PHI at the primary node (I,J).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real ( kind = 8 ) HETA(NI,NJ).
!    HETA is the length in the (X,Y) coordinate system, of the cell
!    in the ETA direction.  This is the actual length of the line
!    connecting the mid-side nodes that separate primary node (I,J)
!    from primary nodes (I,J-1) and (I,J+1).
!
!    Input, real ( kind = 8 ) HKSI(NI,NJ).
!    HKSI is the length in the (X,Y) coordinate system of the cell
!    in the KSI direction.  This is the actual length of the line
!    connecting the mid-side nodes that separate primary node (I,J)
!    from primary nodes (I-1,J) and (I+1,J).
!
!    Input, integer ( kind = 4 ) I, J.
!    I and J are the row and column indices of the
!    primary node at which the gradient of PHI is desired.
!
!    Input, real ( kind = 8 ) PHI(NI,NJ).
!    PHI is the physical quantity whose value is given at each of
!    the primary nodes.
!
!    Output, real ( kind = 8 ) DPHIDX, DPHIDY.
!    DPHIDX and DPHIDY are the computed approximations to the values
!    of the partial derivatives d PHI/d X and d PHI/d Y.
!
!    Input, real ( kind = 8 ) XC(NI,NJ).
!    XC contains the X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Input, real ( kind = 8 ) YC(NI,NJ).
!    YC contains the Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dpdeta
  real ( kind = 8 ) dpdksi
  real ( kind = 8 ) dphidx
  real ( kind = 8 ) dphidy
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdksi
  real ( kind = 8 ) dksidx
  real ( kind = 8 ) dksidy
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydksi
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) phi(ni,nj)
  real ( kind = 8 ) phie
  real ( kind = 8 ) phin
  real ( kind = 8 ) phis
  real ( kind = 8 ) phiw
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) volume
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xe
  real ( kind = 8 ) xn
  real ( kind = 8 ) xs
  real ( kind = 8 ) xw
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ye
  real ( kind = 8 ) yn
  real ( kind = 8 ) ys
  real ( kind = 8 ) yw
!
!  We are given the value of some quantity PHI at the primary nodes.
!
!  We want to estimate the gradients dPHI/dKSI and dPHI/dETA at
!  those primary nodes.
!
!  To do so, we must first estimate the value of PHI at the midpoint nodes,
!  which lie between the given primary node and each of its four neighbors,
!  above, below, to the right and to the left.
!
!  Our task is made a little easier because of the fact that, in the
!  KSI and ETA coordinate system, primary nodes are separated by
!  exactly one unit, (and so are corner nodes, and midside nodes).
!
!  This will make it easy to estimate PHI as a function of KSI and
!  ETA, and hence to get dPHI/dKSI and dPHI/dETA by using a difference
!  of two values.
!
!  However, since we really want to compute dPHI/dX and dPHI/dY,
!  we must make some adjustments to this problem:
!
!    1) the interpolated value of PHI at the midside node must not be
!       the average of the two nearby primary node values, but the WEIGHTED
!       average.  The weights come from the actual X, Y distances between
!       the primary nodes and the midside node in the "real world".
!
!    2) in the X, Y coordinate system, some nodes are separated by
!       zero distance.  This is just to help us handle boundary conditions,
!       but we have to be careful to avoid getting a zero divisor in
!       our weights.
!
!     3) we start out by calculating dPHI/dKSI and dPHI/dETA.  To get
!        dPHI/dX and dPHI/dY, we must compute the jacobian of the
!        transformation that maps KSI and ETA to X and Y.  Then we
!        have to compute the INVERSE of this 2 by 2 matrix, to get
!        the terms dKSI/dX, dKSI/dY, dETA/dX and dETA/dY.  Then we
!        can use the chain rule to calculate terms like:
!
!          dPHI/dX  =  dPHI/dKSI * dKSI/dX + dPHI/dETA * dETA/dX
!
!
!  Now let us look at our primary node in the simple (KSI,ETA)
!  coordinate system, where:
!
!    the (I,J) primary node has (KSI,ETA) coordinates (I,J).
!
!    the midside nodes that neighbor primary node (I,J) have
!    (KSI,ETA) coordinates as follows:
!
!      North: (I, J+1/2)
!      South: (I, J-1/2)
!      East:  (I+1/2, J)
!      West:  (I-1/2, J)
!
!    the corner nodes that form the corners of the cell containing primary
!    node (I,J) have (KSI,ETA) coordinates as follows:
!
!      Northwest: (I-1/2, J+1/2)
!      Northeast: (I+1/2, J+1/2)
!      Southwest: (I-1/2, J-1/2)
!      Southeast: (I+1/2, J-1/2)
!
!    However, if we want to get the X, Y coordinates of a corner node,
!    we can't use these (KSI,ETA) subscripts as indices, because FORTRAN
!    only accepts whole number indices.  So the appropriate indices to
!    use to access the XC, YC arrays are:
!
!      Northwest: (I,   J+1)
!      Northeast: (I+1, J+1)
!      Southwest: (I,   J)
!      Southeast: (I+1, J)
!
!  In the following diagram, we are interested primary node (I,J).
!
!  We show the cell containing (I,J), the neighbor four primary nodes,
!  the corner nodes, and the midside nodes.
!
!
!   ^                   |         (I, J+1)          |
!   |                   |                           |
!   |                   |                           |
!   |   ----------[I-1/2, J+1/2]--{I, J+1/2}--[I+1/2, J+1/2]----
!   |                   |                           |
!   E                   |                           |
!   T   (I-1, J)  {I-1/2, J]      (I, J)      {I+1/2, J}     (I+1/2, J)
!   A                   |                           |
!   |                   |                           |
!   |   ----------[I-1/2, J-1/2]--{I, J-1/2}--[I+1/2, J-1/2]------
!   |                   |                           |
!   |                   |                           |
!   |                   |         (I, J-1)          |
!   |
!   |
!   +-------------------KSI ( = I for primary nodes) ------------->
!
!
!  Make a weighted average of PHI at the east and at the center,
!  to estimate PHIE, the value of PHI at the east midside node:
!
  if ( hksi(i,j)+hksi(i+1,j) /= 0.0D+00 ) then
    t1 = hksi(i,j)/(hksi(i,j)+hksi(i+1,j))
    t2 = hksi(i+1,j)/(hksi(i,j)+hksi(i+1,j))
  else
    t1 = 0.5D+00
    t2 = 0.5D+00
  end if

  phie = t1*phi(i+1,j)+t2*phi(i,j)
!
!  Make a weighted average of PHI at the north and at the center,
!  to estimate PHIN, the value of PHI at the north midside node.
!
  if ( heta(i,j)+heta(i,j+1) /= 0.0D+00 ) then
    t1 = heta(i,j)/(heta(i,j)+heta(i,j+1))
    t2 = heta(i,j+1)/(heta(i,j)+heta(i,j+1))
  else
    t1 = 0.5D+00
    t2 = 0.5D+00
  end if

  phin = t1*phi(i,j+1)+t2*phi(i,j)
!
!  Make a weighted average of PHI at the west and at the center,
!  to estimate PHIW, the value of PHI at the west midside node.
!
  if ( hksi(i,j)+hksi(i-1,j) /= 0.0D+00 ) then
    t1 = hksi(i,j)/(hksi(i,j)+hksi(i-1,j))
    t2 = hksi(i-1,j)/(hksi(i,j)+hksi(i-1,j))
  else
    t1 = 0.5D+00
    t2 = 0.5D+00
  end if

  phiw = t1*phi(i-1,j)+t2*phi(i,j)
!
!  Make a weighted average of PHI at the south and at the center,
!  to estimate PHIS, the value of PHI at the south midside node.
!
  if ( heta(i,j)+heta(i,j-1) /= 0.0D+00 ) then
    t1 = heta(i,j)/(heta(i,j)+heta(i,j-1))
    t2 = heta(i,j-1)/(heta(i,j)+heta(i,j-1))
  else
    t1 = 0.5D+00
    t2 = 0.5D+00
  end if

  phis = t1*phi(i,j-1)+t2*phi(i,j)
!
!  Now subtract opposing values, to estimate dPHI/dKSI and dPHI/dETA
!  at the primary node.  We are assuming that DELTA KSI  =  DELTA ETA = 1.
!
  dpdksi = phie-phiw
  dpdeta = phin-phis
!
!  Now we have to convert our results, which are in terms of KSI and
!  ETA, into the X and Y coordinate system.
!
!  Using the X, Y coordinates of the "corner nodes", compute
!  the X, Y coordinates of the midside nodes of the control volume
!  interfaces that surround the primary node with coordinates (I,J).
!
!    ^
!    |     [I-1/2, J+1/2]-----{I, J+1/2}----[I+1/2, J+1/2]
!    E           |                                |
!    T     {I-1/2, J}         (I, J)        {I+1/2, J}
!    A           |                                |
!    |     [I-1/2, J-1/2]-----{I, J-1/2}----[I+1/2, J]
!    |
!    +----------------------XSI-------------------------->
!
  xe = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j))
  ye = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j))

  xw = 0.5D+00*(xc(i,j+1)+xc(i,j))
  yw = 0.5D+00*(yc(i,j+1)+yc(i,j))

  xn = 0.5D+00*(xc(i,j+1)+xc(i+1,j+1))
  yn = 0.5D+00*(yc(i,j+1)+yc(i+1,j+1))

  xs = 0.5D+00*(xc(i,j)+xc(i+1,j))
  ys = 0.5D+00*(yc(i,j)+yc(i+1,j))
!
!  From the X and Y coordinates of the midside nodes,
!  estimate dX/dKSI, dX/dETA, dY/dKSI and dY/dETA at the primary node,
!  again assuming a spacing delta KSI or delta ETA  =  1
!  for opposing midside nodes.
!
  dxdksi = xe-xw
  dxdeta = xn-xs
  dydksi = ye-yw
  dydeta = yn-ys
!
!  Compute the determinant of the Jacobian (XSI,ETA)-->(X,Y).
!
!  J(XSI,ETA)  =  ( dX/dKSI  dX/dETA)
!               ( dY/dKSI  dY/dETA)
!
  volume = dxdksi*dydeta-dxdeta*dydksi
!
!  If the determinant is zero, then the control volume has zero
!  volume, which should never happen.
!
  if ( volume == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GRADNT - Fatal error!'
    write ( *, '(a)' ) '  The control volume has zero volume!'
    stop
  end if
!
!  Now compute the elements of the inverse Jacobian.
!
!  J_Inverse(XSI,ETA)  =  ( dKSI/dX  dKSI/dY)
!                       ( dETA/dX  dETA/dY)
!
!  using the simple fact that, for a 2 by 2 matrix,
!
!  A_Inverse  =    (a22  -a21)  / determinant(A)
!                (-a12  a11)
!
  dksidx = dydeta/volume
  dksidy = -dxdeta/volume
  detadx = -dydksi/volume
  detady = dxdksi/volume
!
!  Now we can finally compute dPHI/dX and dPHI/dY using the chain rule.
!
  dphidx  =  dpdksi*dksidx + dpdeta*detadx
  dphidy  =  dpdeta*detady + dpdksi*dksidy

  return
end
subroutine hello ( liwork, lwork )

!*****************************************************************************80
!
!! HELLO prints a brief message identifying the program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ) liwork
  integer ( kind = 4 ) lwork

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CRYSTAL_QED:'
  write ( *, '(a)' ) '  Version of 29 October 1996'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  CRYSTAL_QED seeks a set of parameters that minimize a '
  write ( *, '(a)' ) '  cost functional associated with a crystallization '
  write ( *, '(a)' ) '  problem.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  LIWORK, integer workspace size  =  ',liwork
  write ( *, '(a,i6)' ) '  LWORK, real workspace size  =      ',lwork

  return
end
subroutine inidat(ae1,ae2,ak1,ak2,area,b,b1jbl,b2jbl,b3jbl, &
    birad,bo,cappa,cd,ce1,ce2,cfo,cinc,cmu,cost,cvn,delt,dtm, &
    epsad,epsil,ewall,f,fcsl,fks,fksl,fma,fn,fnu,fo,fr,frsl, &
    gam,gamt,grash,hamag,heta,hf,hksi,icost,icrys,inturb,iplot, &
    ipref,iprint,jcrys,jpref,l0,l1,last,lastt,lblk,lortho, &
    lsolve,m0,m1,maxbot,mode,nbot,ndt,ni,nj,nk,np,npc,ns, &
    nsolve,ntimes,orth,pr,ra,rdtm, &
    re,recb,rect,relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk, &
    sigma,sigt,smooth,stel,stes,tal,tanca,tanca2,tas,tend,tf,tinit, &
    title,tnow,tw,vave,vol,vort,x,xbot,xc,xlen,y,ybot,yc,ylen)

!*****************************************************************************80
!
!! INIDAT sets the initial values of certain data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Output, real ( kind = 8 ) XC(NI,NJ).
!    The X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Output, real ( kind = 8 ) YC(NI,NJ).
!    The Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ) maxbot
  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj
  integer ( kind = 4 ) nk
  integer ( kind = 4 ) ns

  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) b
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) b3jbl(ni)
  real ( kind = 8 ) birad
  real ( kind = 8 ) bo
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cfo
  real ( kind = 8 ) cinc
  real ( kind = 8 ) cl
  real ( kind = 8 ) cmu
  real ( kind = 8 ) cost
  real ( kind = 8 ) cs
  real ( kind = 8 ) cvn
  real ( kind = 8 ) delt
  real ( kind = 8 ) dtm
  real ( kind = 8 ) epsad
  real ( kind = 8 ) epsil
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fkl
  real ( kind = 8 ) fks
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fnu
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) grav
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hf
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icost
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) inturb
  integer ( kind = 4 ) iplot
  integer ( kind = 4 ) ipref
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) jpref
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lastt
  logical lblk(ns)
  logical lortho
  logical lsolve(ns)
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nbot
  integer ( kind = 4 ) ndt
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npc
  integer ( kind = 4 ) nsolve(ns)
  integer ( kind = 4 ) ntimes(ns)
  real ( kind = 8 ) orth
  real ( kind = 8 ) pr
  real ( kind = 8 ) ra
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) rhol
  real ( kind = 8 ) rhos
  real ( kind = 8 ) rpr
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigma
  real ( kind = 8 ) sigt
  real ( kind = 8 ) smooth
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) tal
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) tas
  real ( kind = 8 ) tend
  real ( kind = 8 ) tf
  real ( kind = 8 ) tin
  real ( kind = 8 ) tinit
  character ( len = 25 ) title(ns)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) tw
  real ( kind = 8 ) vave
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) vort(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xbot(maxbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ybot(maxbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ylen
!
!  Set quantities whose values, initial values, or dummy values
!  do not depend on other quantities.
!
  do i = 1,ni
    do j = 1,nj
      ae1(i,j) = 0.0D+00
      ae2(i,j) = 0.0D+00
      ak1(i,j) = 0.0D+00
      ak2(i,j) = 0.0D+00
    end do
  end do

  do i = 1,ni
    do j = 1,nj
      area(i,j) = 0.0D+00
    end do
  end do

  b = 0.1778D+00

  do i = 1,nj
    b1jbl(i) = 0.0D+00
    b2jbl(i) = 0.0D+00
  end do

  do i = 1,ni
    b3jbl(i) = 0.0D+00
  end do

  cappa = 0.4187D+00
  cd = 1.0D+00
  ce1 = 1.44D+00
  ce2 = 1.92D+00
  cinc = 0.0D+00
  cl = 1000.0D+00
  cmu = 0.09D+00
  cost = 0.0D+00
  cs = 1000.0D+00
  cvn = 200000.0D+00
  dtm = 100.0D+00
  epsad = 0.0001D+00
  epsil = 0.0001D+00
  ewall = 9.793D+00

  do i = 1,ni
    do j = 1,nj
      do k = 1,ns
        f(i,j,k) = 0.0D+00
      end do
    end do
  end do

  fkl = 64.0D+00
  fks = 22.0D+00
  fma = -1000.0D+00

  do i = 1,ni
    do j = 1,nj
      do k = 1,ns
        fn(i,j,k) = 0.0D+00
      end do
    end do
  end do

  fnu = 3.0D-07

  do i = 1,ni
    do j = 1,nj
      do k = 1,ns
        fo(i,j,k) = 0.0D+00
      end do
    end do
  end do

  fr = 1.0D-10

  do i = 1,ni
    do j = 1,nj
      gam(i,j) = 0.0D+00
    end do
  end do

  grash = 10000000.0D+00
  grav = 9.81D+00
  hamag = 0.0D+00
  heta(1:ni,1:nj) = 0.0D+00
  hf = 1800000.0D+00

  do i = 1,ni
    do j = 1,nj
      hksi(i,j) = 0.0D+00
    end do
  end do
!
!  ICOST = 1, U**2+V**2
!        2, (T-TF)**2
!        3, (VAVE-SQRT(U**2+V**2))
!
  icost = 1

  icrys = 42
  inturb = 0
!
!  IPLOT = 0, do not create a plot file.
!        1, create plot file.
!
  iplot = 1
  ipref = 11
!
!  IPRINT = 0, don't print very much.
!         1, print intermediate information.
!
  iprint = 0

  jcrys = 22
  jpref = 30
  l0 = 62
  l1 = 62
!
!  Interstep iteration number.
!  Use LAST = 40 for more accurate results.
!  Use LAST = 10 for quicker results.
!
  last = 40
!
!  Set the number of timesteps.
!
!  For tests, set LASTT = 2.
!  For a real run, try LASTT = 200.
!
  lastt = 2

  do i = 1,4
    lblk(i) = .false.
  end do
  do i = 5,10
    lblk(i) = .true.
  end do

  lortho = .false.

  do i = 1,ns
    lsolve(i) = .false.
  end do

  m0 = 42
  m1 = 42
  mode = 1
  ndt = 0
  np = 3
  npc = 4
!
!  Number of linear iterations.
!
  nsolve(1) = 3
  nsolve(2) = 3
  nsolve(3) = 1
  nsolve(4) = 1
  nsolve(5) = 3
  nsolve(6) = 3
  nsolve(7) = 3
  nsolve(8) = 3
  nsolve(9) = 1
  nsolve(10) = 1
!
!  Number of nonlinear iterations.
!
  ntimes(1) = 5
  ntimes(2) = 5
  ntimes(3) = 1
  ntimes(4) = 1
  ntimes(5) = 3
  ntimes(6) = 3
  ntimes(7) = 3
  ntimes(8) = 3
  ntimes(9) = 1
  ntimes(10) = 1

  orth = 50000.0D+00
  pr = 0.015D+00
  rdtm = 0.0D+00
  re = 1.0D+00
  recb = 0.0D+00
  rect = 0.0D+00

  relax(1) = 0.3D+00
  relax(2) = 0.3D+00
  relax(3) = 0.8D+00
  relax(4) = 1.0D+00
  relax(5) = 0.7D+00
  relax(6) = 0.2D+00
  relax(7) = 0.5D+00
  relax(8) = 0.5D+00
  relax(9) = 1.0D+00
  relax(10) = 1.0D+00
  relax(11) = 1.0D+00
  relax(12) = 0.6D+00
  relax(13) = 1.0D+00
  relax(14) = 1.0D+00

  rhocon = 1.0D+00
  rhol = 2490.0D+00
  rhos = 2490.0D+00

  do i = 1,ni
    do j = 1,nj
      rueta(i,j) = 0.0D+00
      ruksi(i,j) = 0.0D+00
    end do
  end do

  sige = 1.3D+00
  sigk = 1.0D+00
  sigma = 0.72D+00
  sigt = 0.9D+00
  smooth = 0.0001D+00
  tal = 1523.0D+00
  tanca = 1.0D+00
  tanca2 = 1.0D+00
  tas = 1523.0D+00
  tf = 1683.0D+00
  tin = 1523.0D+00
  tinit = 0.0D+00

  title(1) = 'U velocity'
  title(2) = 'V velocity'
  title(3) = 'Pressure'
  title(4) = 'Corrected pressure'
  title(5) = 'Temperature'
  title(6) = 'Rotational velocity'
  title(7) = 'Turbulent dissipation'
  title(8) = 'Turbulent energy'
  title(9) = 'Magnetic stream function'
  title(10) = 'Stream function'

  tw = 1713.0D+00
  vave = 1850.0D+00

  do i = 1,ni
    do j = 1,nj
      vol(i,j) = 0.0D+00
    end do
  end do

  do i = 1,ni
    do j = 1,nj
      vort(i,j) = 0.0D+00
    end do
  end do

  do i = 1,ni
    do j = 1,nj
      x(i,j) = 0.0D+00
    end do
  end do

  do i = 1,ni
    do j = 1,nj
      xc(i,j) = 0.0D+00
    end do
  end do

  xlen = 1.2

  do i = 1,ni
    do j = 1,nj
      y(i,j) = 0.0D+00
      yc(i,j) = 0.0D+00
    end do
  end do

  ylen = 1.0
!
!  Set things that depend on things.
!
  birad = 0.25D+00 * 5.67D+00 * 1.5D+00 * 0.0254D+00 /fks/(tw-tf)/re
  bo = (rhol*grav*b**2)/sigma
  cfo = fnu/(b*b)

  if ( inturb == 1 ) then
    do i = 1,ni
      do j = 1,nj
        f(i,j,7) = 0.005D+00
      end do
    end do
  end if

  fcsl = cs/cl
  fksl = fks/fkl
  frsl = rhos/rhol
  ra = grash*pr
  rho(1:ni,1:nj) = rhocon
  rpr = rhocon/(pr*re)
  stel = cl*(tw-tf)/hf
  stes = cs*(tf-tin)/hf
  tend = tinit+dtm*lastt
  tnow = tinit
!
!  Set things that depend on things that depend on things.
!
  delt = cfo*dtm

  if ( inturb == 1 ) then
    do i = 1,ni
      do j = 1,nj
        f(i,j,8) = f(i,j,7)**1.5 / 0.006D+00
        gamt(i,j) = 0.006D+00 *cmu*rho(i,j)*f(i,j,7)**0.5
      end do
    end do
  end if
!
!  Set the crucible shape.
!
  xbot(1) = 0.0D+00
  ybot(1) = 0.0D+00

  xbot(2) = 0.015D+00
  ybot(2) = 0.3D+00

  xbot(3) = 0.03D+00
  ybot(3) = 0.5D+00

  xbot(4) = 0.046D+00
  ybot(4) = 0.6D+00

  xbot(5) = 0.07D+00
  ybot(5) = 0.7D+00

  xbot(6) = 0.10D+00
  ybot(6) = 0.8D+00

  xbot(7) = 0.14D+00
  ybot(7) = 0.9D+00

  xbot(8) = 0.18D+00
  ybot(8) = 0.95D+00

  xbot(9) = 0.25D+00
  ybot(9) = 1.0D+00

  nbot = 9

  return
end
subroutine inigrd ( icrys, jcrys, l0, m0, nbot, xbot, xc, xlen, ybot, &
  yc, ylen )
!
!*****************************************************************************80
!
!! INIGRD makes an initial assignment of the grid points XC, YC.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Output, real ( kind = 8 ) XC(NI,NJ).
!    The X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Output, real ( kind = 8 ) YC(NI,NJ).
!    The Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ) nbot
  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  real ( kind = 8 ) free
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) licrys
  integer ( kind = 4 ) ll1
  integer ( kind = 4 ) llst1
  integer ( kind = 4 ) m0
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) wave
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) ylen

  ll1 = icrys/2+1
  llst1 = jcrys/2+1
  licrys = (jcrys+m0)/2+1

  free = 0.5D+00*xlen
!
!  Set the Y coordinates
!
  do i = 2, l0

    if ( i < ll1 ) then
      t2 = free*(0.5D+00*(dble(i-2)/dble(ll1-2))**1.5)
    else if ( i < icrys ) then
      t2 = free*(1.0D+00-0.5D+00*(dble(icrys-i)/dble(icrys-ll1))**1.5)
    else
      t2 = free+(xlen-free)*(dble(i-icrys)/dble(l0-icrys))**1.2
    end if

    do j = 2, m0

      wave = 0.4D+00 + 0.02D+00 * sin ( t2 * 20.0D+00 * 3.14159D+00 )

      if ( j < llst1 ) then
        t1 = wave/2.0D+00*(dble(j-2)/dble(llst1-2))**1.5
      else if ( j < jcrys ) then
        t1 = wave-wave/2.0D+00*(dble(jcrys-j)/dble(jcrys-llst1))**1.5
      else if ( j < licrys ) then
        t1 = 0.4D+00 + 0.3D+00*(dble(j-jcrys)/dble(licrys-jcrys))**1.5
      else
        t1 = 1.0D+00 - 0.3D+00*(dble(m0-j)/dble(m0-licrys))**1.5
      end if

      yc(i,j) = t1*ylen

    end do
  end do
!
!  Set the XC coordinates along the bottom boundary.
!
  xc(2,2) = xbot(1)

  do j = 3, m0-1

    call cubic ( nbot, yc(2,j), ybot, xc(2,j), xbot )

  end do

  xc(2,m0) = xbot(nbot)
!
!  Now assign the XC coordinates of the other nodes.
!  Here, we explicitly assume that the free surface and crystal
!  boundaries occur at a coordinate value of 0.6.
!
  do j = 2, m0
    do i = 3, l0

      if ( i < ll1 ) then

        t1 = 0.5D+00 * (dble(i-2)/dble(ll1-2))**1.5
        xc(i,j) = (1.0D+00 - t1)*xc(2,j)+t1*free

      else if ( i <= icrys ) then

        t1 = 1.0D+00 - 0.5D+00 * (dble(icrys-i)/dble(icrys-ll1))**1.5
        xc(i,j) = ( 1.0D+00 - t1)*xc(2,j)+t1*free

      else

        t1 = (dble(i-icrys)/dble(l0-icrys))**1.2

        if ( j < jcrys ) then
          xc(i,j) = xc(icrys,j)+t1*((xlen-free)+0.5D+00*yc(l0,jcrys)**2 &
            -0.5D+00 * yc(l0,j)**2)
        else
          xc(i,j) = xc(icrys,j)+t1*(xlen-free)
        end if

      end if

    end do
  end do
!
!  Copy values to dummy nodes with I = 1 or J=1.
!
  xc(1,1) = xc(2,2)
  yc(1,1) = yc(2,2)

  do j = 2, m0
    xc(1,j) = xc(2,j)
    yc(1,j) = yc(2,j)
  end do

  do i = 2, l0
    xc(i,1) = xc(i,2)
    yc(i,1) = yc(i,2)
  end do

  return
end
subroutine movgrd ( b1jbl, b2jbl, bo, delt, fksl, fr, frsl, icrys, &
  iprint, jcrys, l0, m0, p, pr, re, stel, tanca, tanca2, xc, xlen, yc )

!*****************************************************************************80
!
!! MOVGRD calculates the new position of the interface and free surface.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input/output, real ( kind = 8 ) XC(NI,NJ).
!    The X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Input/output, real ( kind = 8 ) YC(NI,NJ).
!    The Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  real ( kind = 8 ) ac
  real ( kind = 8 ) as
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) bo
  real ( kind = 8 ) coeff
  real ( kind = 8 ) coff
  real ( kind = 8 ) d2hdy2(nj)
  real ( kind = 8 ) delt
  real ( kind = 8 ) dhdy(nj)
  real ( kind = 8 ) dhdym(nj)
  real ( kind = 8 ) dhdyp(nj)
  real ( kind = 8 ) dlen
  real ( kind = 8 ) dx1
  real ( kind = 8 ) dxm1
  real ( kind = 8 ) dy1
  real ( kind = 8 ) dym1
  real ( kind = 8 ) fksl
  real ( kind = 8 ) flow(ni)
  real ( kind = 8 ) flomax
  real ( kind = 8 ) flomin
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) high(ni)
  real ( kind = 8 ) hilev
  real ( kind = 8 ) himax
  real ( kind = 8 ) himin
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) interl
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) m0
  real ( kind = 8 ) omega
  real ( kind = 8 ) omega2
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) pcsurf
  real ( kind = 8 ) pr
  real ( kind = 8 ) pres(ni)
  real ( kind = 8 ) pullv
  real ( kind = 8 ) re
  real ( kind = 8 ) rr(nj)
  real ( kind = 8 ) stel
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) xb(ni)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xlen
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) zb1(ni)

  do j = 2,m0
    zb1(j) = xc(icrys,j)
  end do

  do k = 1,5

    do j = jcrys+1,m0-1
      dhdy(j) = (xc(icrys,j+1)-xc(icrys,j))/(yc(icrys,j+1)-yc(icrys,j))
      dhdym(j) = (xc(icrys,j)-xc(icrys,j-1))/(yc(icrys,j)-yc(icrys,j-1))
      dhdyp(j) = (xc(icrys,j+1)-xc(icrys,j))/(yc(icrys,j+1)-yc(icrys,j))
    end do

    do j = jcrys+1,m0-1
      d2hdy2(j) = (dhdyp(j)-dhdym(j))/(yc(icrys,j+1)-yc(icrys,j))
      rr(j) = d2hdy2(j)/(1.0D+00 + (dhdy(j)*dhdy(j)))**1.5 &
        +dhdy(j)/yc(icrys,j)/(1.0D+00 + (dhdy(j)*dhdy(j)))**0.5
    end do

    pres(jcrys+1) = 0.0D+00

    do j = jcrys+1,m0-2
      pres(j+1) = pres(j)+(p(icrys-1,j+1)-p(icrys-1,j))/yc(icrys-1,j+1)
    end do

    himax = 0.0D+00
    himin = 0.0D+00
    interl = (jcrys+m0)/2

    do j = jcrys+1,m0-1
      pcsurf = pres(j)-pres(interl)
      high(j) = fr*pcsurf+rr(j)/bo
    end do

    hilev = high(interl)

    omega = 0.0D+00

    do j = jcrys+1,m0-1
      high(j) = omega*((high(j)-hilev)-(xc(icrys,j)-xc(icrys,interl)))
      himax = max(himax,high(j))
      himin = min(himin,high(j))
      xc(icrys,j) = xc(icrys,j)+high(j)
    end do

    if ( omega /= 0.0 ) then
      xc(icrys,jcrys) = xc(icrys,jcrys+1) &
        +tanca*(yc(icrys,jcrys+1)-yc(icrys,jcrys))
      xc(icrys,m0) = xc(icrys,m0-1)+tanca2*(yc(icrys,m0)-yc(icrys,m0-1))
    end if

    flomax = 0.0D+00
    flomin = 0.0D+00

    if ( iprint > 0 ) then
      if ( k == 5 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MOVGRD:'
        write ( *, '(a,g14.6)' ) '  Himin = ',himin
        write ( *, '(a,g14.6)' ) '  Himax = ',himax
        write ( *, '(a,g14.6)' ) '  Omega = ',omega
        write ( *, '(a)' ) ' '
        write(*,'(7f10.4)')high(jcrys+1),high(jcrys+2),high(m0-10), &
          high(m0-8),high(m0-4),high(m0-2),high(m0-1)
        write(*,'(7f9.4)')xc(icrys,jcrys),xc(icrys,jcrys+1), &
          xc(icrys,jcrys+4),xc(icrys,m0-8),xc(icrys,m0-4), &
          xc(icrys,m0-1),xc(icrys,m0)
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MOVGRD:'
        write ( *, '(a)' ) '  Free surface movement'
      end if
    end if

  end do

  omega2 = 0.0D+00

  do j = 2,jcrys-1

    dxm1 = xc(icrys,j+1)-xc(icrys,j)
    dym1 = yc(icrys,j+1)-yc(icrys,j)
    dlen = sqrt(dxm1*dxm1+dym1*dym1)
    as = -dxm1/dlen
    ac = dym1/dlen
    dx1 = delt*frsl*(b1jbl(j)-fksl*b2jbl(j))*stel/(re*pr)*ac
    dy1 = delt*frsl*(b1jbl(j)-fksl*b2jbl(j))*stel/(re*pr)*as

    if ( j == jcrys-1 ) then
      flow(j) = dx1
    else
      flow(j) = dx1+dy1**2/dx1
    end if

  end do
!
!  What is PULLV?
!
  pullv = flow(jcrys-1)

  do j = 3,jcrys-1
    flomax = max(flomax,(flow(j)-pullv))
    flomin = min(flomin,(flow(j)-pullv))
  end do

  omega2 = min(omega2, 0.003D+00/(flomax+1.0D-06))
  omega2 = min(omega2, 0.003D+00/(-flomin+1.0D-06))

  do j = 3,jcrys-1
    flow(j) = omega2*(flow(j)-pullv)+(xc(icrys,jcrys)-zb1(jcrys))
    xb(j) = xc(icrys,j)+flow(j)
    xc(icrys,j) = xb(j)
  end do

  xc(icrys,2) = xc(icrys,3)

  if ( iprint > 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOVGRD:'
    write ( *, '(a,g14.6)' ) '  Flomin  =  ',flomin
    write ( *, '(a,g14.6)' ) '  Flomax  =  ',flomax
    write ( *, '(a,g14.6)' ) '  Omega  =   ',omega2
    write ( *, '(a)' ) ' '
    write(*,'(2x,7f9.4)')xb(jcrys-1),xb(jcrys-2),xb(8),xb(6),xb(4),xb(3),xb(2)
    write(*,'(2x,7f9.4)')b1jbl(jcrys-1),-fksl*b2jbl(jcrys-1), &
      b1jbl(jcrys-2),-fksl*b2jbl(jcrys-2),b1jbl(jcrys-5),-fksl*b2jbl(jcrys-5)
    write(*,'(2x,7f9.4)') b1jbl(jcrys-1)-fksl*b2jbl(jcrys-1), &
      b1jbl(jcrys-2)-fksl*b2jbl(jcrys-2),b1jbl(jcrys-5)-fksl*b2jbl(jcrys-5), &
      b1jbl(jcrys-10)-fksl*b2jbl(jcrys-10),b1jbl(jcrys-15)-fksl*b2jbl(jcrys-15), &
      b1jbl(jcrys-17)-fksl*b2jbl(jcrys-17),b1jbl(jcrys-19)-fksl*b2jbl(jcrys-19)
    write(*,'(2x,7f9.4)')flow(21),flow(19),flow(17),flow(15), &
      flow(13),flow(11),flow(10)
    write(*,'(2x,7f9.4)')flow(9),flow(8),flow(7),flow(6),flow(5),flow(4),flow(3)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'MOVGRD:'
    write ( *, '(a)' ) '  Solid-liquid interface movement'
    write ( *, '(a)' ) ' '
  end if

  do j = 2,m0

    coeff = (xc(2,j)-xc(icrys,j))/(xc(2,j)-zb1(j))
    coff = (xlen-xc(icrys,j))/(xlen-zb1(j))

    do i = 2,icrys-1
      xc(i,j) = xc(2,j)-(xc(2,j)-xc(i,j))*coeff
    end do

    do i = icrys+1,l0
      xc(i,j) = xlen-(xlen-xc(i,j))*coff
      if ( j <= jcrys .and. i == l0 ) then
        xc(l0,j) = xlen+0.5D+00*yc(l0,jcrys)**2 - 0.5D+00*yc(l0,j)**2
      end if
    end do

  end do

  return
end
subroutine output ( cfo, iter, izone, res, smax, ssum, t, tnow, u, w )

!*****************************************************************************80
!
!! OUTPUT prints some sample data from the ongoing solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: ns = 10

  real ( kind = 8 ) cfo
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) izone
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) smax
  real ( kind = 8 ) ssum
  real ( kind = 8 ) t(ni,nj)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) w(ni,nj)
!
!  Solid zone output.
!
  if ( izone == 1 ) then

    if ( iter == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OUTPUT:'
      write ( *, '(a,g14.6)' ) '  Current time  =  ',tnow
      write ( *, '(a,g14.6,a)' ) '   =  ',tnow*cfo,' seconds.'
      write ( *, '(a)' ) '  Solid zone results.'
      write ( *, '(a)' ) '  Iter    Res(5)      T(50,21)    T(50,22)'
      write ( *, '(a)' ) ' '
    end if

    write(*,'(i4,2x,10g12.4)')iter,res(5),t(50,21),t(50,22)
!
!  Liquid zone output.
!
  else if ( izone == 2 ) then

    if ( iter == 1 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OUTPUT:'
      write ( *, '(a)' ) '  Liquid zone results.'
      write ( *, '(a)' ) '  SMAX is the maximum local mass imbalance.'
      write ( *, '(a)' ) '  SSUM is the total mass imbalance.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' Iter    SMAX        SSUM        U(5,5)      W(5,5)      T(5,5)'
      write ( *, '(a)' ) ' '
    end if

    write(*,'(i4,2x,10g12.4)')iter,smax,ssum,u(5,5),w(5,5),t(5,5)

  end if

  return
end
subroutine pmod ( ipref, jcrys, jpref, l1, m1, p )

!*****************************************************************************80
!
!! PMOD interpolates the pressure at boundary corners.
!
!  Discussion:
!
!    PMOD carries out some simple operations to extend the value of the
!    pressure to primary nodes at the corners, and to normalize the
!    pressure by subtracting off its value at the reference point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipref
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) jpref
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) pref
!
!  Extrapolate to get pressures on the boundary corners.
!
  p(1,1) = p(2,1)+p(1,2)-p(2,2)
  p(l1,1) = p(l1-1,1)+p(l1,2)-p(l1-1,2)
  p(1,jcrys) = p(2,jcrys)+p(1,jcrys-1)-p(2,jcrys-1)
  p(l1,jcrys) = p(l1-1,jcrys)+p(l1,jcrys-1)-p(l1-1,jcrys-1)
!
!  Subtract off the reference pressure.
!
  pref = p(ipref,jpref)

  do j = 1,m1
    do i = 1,l1
      p(i,j) = p(i,j)-pref
    end do
  end do

  return
end
subroutine prdat(b,birad,bo,cappa,cfo,cvn,delt,dtm,fcsl,fksl, &
  fma,fnu,fr,frsl,grash,icost,icrys,inturb,iprint, &
  jcrys,l0,last,lastt,m0,mode,nbot,ns,nsolve,ntimes,orth, &
  pr,ra,rdtm,recb,rect,rhocon,smooth,stel,stes,tanca,tanca2,tend, &
  tf,tinit,title,tw,vave,xlen,ylen)

!*****************************************************************************80
!
!! PRDAT prints out the initial values of certain data.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ) ns

  real ( kind = 8 ) b
  real ( kind = 8 ) birad
  real ( kind = 8 ) bo
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cfo
  real ( kind = 8 ) cvn
  real ( kind = 8 ) delt
  real ( kind = 8 ) dtm
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fksl
  real ( kind = 8 ) fma
  real ( kind = 8 ) fnu
  real ( kind = 8 ) fr
  real ( kind = 8 ) frsl
  real ( kind = 8 ) grash
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icost
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) inturb
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lastt
  integer ( kind = 4 ) m0
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) nbot
  integer ( kind = 4 ) nsolve(ns)
  integer ( kind = 4 ) ntimes(ns)
  real ( kind = 8 ) orth
  real ( kind = 8 ) pr
  real ( kind = 8 ) ra
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) smooth
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) tanca
  real ( kind = 8 ) tanca2
  real ( kind = 8 ) tend
  real ( kind = 8 ) tf
  real ( kind = 8 ) tinit
  character ( len = 25 ) title(ns)
  real ( kind = 8 ) tw
  real ( kind = 8 ) vave
  real ( kind = 8 ) xlen
  real ( kind = 8 ) ylen

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRDAT:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  B, Crucible radius  =                 ',b
  write ( *, '(a,g14.6)' ) '  BIRAD, thermal coefficient  =         ',birad
  write ( *, '(a,g14.6)' ) '  BO, the Bond number  =                ',bo
  write ( *, '(a,g14.6)' ) '  CAPPA, wall function constant  =      ',cappa
  write ( *, '(a,g14.6)' ) '  CFO, time scale  =                    ',cfo
  write ( *, '(a,g14.6)' ) '  CVN, ADAPT cell volume weight  =      ',cvn
  write ( *, '(a,g14.6)' ) '  DELT, time step in seconds  =         ',delt
  write ( *, '(a,g14.6)' ) '  DTM, dimensionless time step  =       ',dtm
  write ( *, '(a,g14.6)' ) '  FCSL  =  CS/CL =                      ',fcsl
  write ( *, '(a,g14.6)' ) '  FKSL  =  FKS/FKL =                    ',fksl
  write ( *, '(a,g14.6)' ) '  FMA, Marangoni number FMA  =	      ',fma
  write ( *, '(a,g14.6)' ) '  FR, Froude number  =		      ',fr
  write ( *, '(a,g14.6)' ) '  FRSL  =  RHOS/RHOL =                  ',frsl
  write ( *, '(a,g14.6)' ) '  GRASH, Grashof number  =              ',grash
  write ( *, '(a,i6)' ) '  ICOST, cost functional  =             ',icost
  write ( *, '(a,i6)' ) '  ICRYS, maximum I of crystal  =        ',icrys
  write ( *, '(a,i6)' ) '  INTURB, turbulence option  =          ',inturb
  write ( *, '(a,i6)' ) '  IPRINT, printing option  =            ',iprint
  write ( *, '(a,i6)' ) '  JCRYS, maximum J of crystal  =        ',jcrys
  write ( *, '(a,i6)' ) '  L0, number of I nodes  =              ',l0
  write ( *, '(a,i6)' ) '  LAST, number of zone iterations on  '
  write ( *, '(a,i6)' ) '    each time step  =                   ',last
  write ( *, '(a,i6)' ) '  LASTT, number of time steps        =  ',lastt
  write ( *, '(a,i6)' ) '  M0, number of J nodes  =              ',m0
  write ( *, '(a,i6)' ) '  MODE, 0 cartesian, 1 axisymmetric  =  ',mode
  write ( *, '(a,i6)' ) '  NBOT, number of boundary points  =    ',nbot
  write ( *, '(a,g14.6)' ) '  ORTH, ADAPT orthogonality weight  =   ',orth
  write ( *, '(a,g14.6)' ) '  PR, Prandtl number PR  =              ',pr
  write ( *, '(a,g14.6)' ) '  RA, Rayleigh number  =                ',ra
  write ( *, '(a,g14.6)' ) '  RBD  =  FMA/RA =                      ',fma/ra
  write ( *, '(a,g14.6)' ) '  RDTM,  =  1/DTM or 0 =                ',rdtm
  write ( *, '(a,g14.6)' ) '  RECB, crucible Reynolds number  =     ',recb
  write ( *, '(a,g14.6)' ) '  RECT, crystal Reynolds number  =      ',rect
  write ( *, '(a,g14.6)' ) '  RHOCON, density constant  =           ',rhocon
  write ( *, '(a,g14.6)' ) '  SMOOTH, ADAPT smoothness weight  =    ',smooth
  write ( *, '(a,g14.6)' ) '  STEL, liquid Stefan number STEL  =    ',stel
  write ( *, '(a,g14.6)' ) '  STES, solid Stefan number STES  =     ',stes
  write ( *, '(a,g14.6)' ) '  TANCA  =                              ',tanca
  write ( *, '(a,g14.6)' ) '  TANCA2  =                             ',tanca2
  write ( *, '(a,g14.6)' ) '  TEND, dimensionless end time  =       ',tend
  write ( *, '(a,g14.6)' ) '  TF, crystal melting temperature  =    ',tf
  write ( *, '(a,g14.6)' ) '  TINIT, dimensionless start time  =    ',tinit
  write ( *, '(a,g14.6)' ) '  TW, the wall temperature  =           ',tw
  write ( *, '(a,g14.6)' ) '  VAVE, desired average velocity  =     ',vave
  write ( *, '(a,g14.6)' ) '  WEBER  =  FR*BO =                     ',fr*bo
  write ( *, '(a,g14.6)' ) '  XLEN  =  problem region length =      ',xlen
  write ( *, '(a,g14.6)' ) '  YLEN  =  problem region height =      ',ylen
  write ( *, '(a,g14.6)' ) '  Characteristic velocity FNU/B  =      ',fnu/b
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Variable                     Nonlinear    Linear'
  write ( *, '(a)' ) '                             Iterations   Iterations'
  write ( *, '(a)' ) ' '
  do i = 1,ns
    write ( *, '(a,i6,i6)' ) title(i),ntimes(i),nsolve(i)
  end do

  return
end
subroutine resid ( aim, aip, ajm, ajp, ap, con, f, l1, m1, nf )

!*****************************************************************************80
!
!! RESID computes the residual of the linear equations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: ns = 10

  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) conmax
  real ( kind = 8 ) diamax
  real ( kind = 8 ) f(ni,nj,ns)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) nf
  real ( kind = 8 ) offmax
  real ( kind = 8 ) ratio
  real ( kind = 8 ) ratmin
  real ( kind = 8 ) res
  real ( kind = 8 ) resmax
  real ( kind = 8 ) sum
  real ( kind = 8 ) xmax

  iprint = 0

  conmax = 0.0D+00
  diamax = 0.0D+00
  offmax = 0.0D+00
  ratmin = 100000.0D+00
  resmax = 0.0D+00
  xmax = 0.0D+00

  imax = 0
  jmax = 0

  do i = 2,l1-1
    do j = 2,m1-1

      res = con(i,j)-ap(i,j)*f(i,j,nf)+aim(i,j)*f(i-1,j,nf) &
        +aip(i,j)*f(i+1,j,nf) &
        +ajm(i,j)*f(i,j-1,nf)+ajp(i,j)*f(i,j+1,nf)

      if ( abs(res) >= resmax ) then
        resmax = abs(res)
        imax = i
        jmax = j
      end if

      if ( abs(con(i,j)) >= conmax ) then
        conmax = abs(con(i,j))
      end if

      if ( abs(f(i,j,nf)) >= xmax ) then
        xmax = abs(f(i,j,nf))
      end if

      if ( abs(ap(i,j)) >= diamax ) then
        diamax = abs(ap(i,j))
      end if

      sum = abs(aim(i,j))+abs(aip(i,j))+abs(ajm(i,j))+abs(ajp(i,j))
      if ( sum >= offmax ) then
        offmax = sum
      end if

      if ( sum /= 0.0 ) then
        ratio = abs(ap(i,j))/sum
        if ( ratio < ratmin ) then
          ratmin = ratio
        end if
      end if

    end do
  end do

  if ( iprint == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RESID'
    write ( *, '(a,i6)' ) '  Maxixum residual for variable ',nf
    write ( *, '(a,g14.6)' ) '  is ',resmax
    write ( *, '(a,2i6)' ) '  at I,J = ',imax,jmax
    i = imax
    j = jmax
    write ( *, '(a)' ) ' '
    write ( *, '(a,2g14.6)' ) '  A(I,J)*X(I,J) =     ',ap(i,j),f(i,j,nf)
    write ( *, '(a,2g14.6)' ) '  A(I-1,J)*X(I-1,J) = ',aim(i,j),f(i-1,j,nf)
    write ( *, '(a,2g14.6)' ) '  A(I+1,J)*X(I+1,J) = ',aip(i,j),f(i+1,j,nf)
    write ( *, '(a,2g14.6)' ) '  A(I,J-1)*X(I,J-1) = ',ajm(i,j),f(i,j-1,nf)
    write ( *, '(a,2g14.6)' ) '  A(I,J+1)*X(I,J+1) = ',ajp(i,j),f(i,j+1,nf)
    write ( *, '(a,g14.6)' ) '  CON(I,J) =          ',con(i,j)
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Norm of variable is ',xmax
    write ( *, '(a,g14.6)' ) '  Norm of RHS is      ',conmax
    write ( *, '(a,g14.6)' ) '  Norm of diagonal is ',diamax
    write ( *, '(a,g14.6)' ) '  Norm of off-diag is ',offmax
    write ( *, '(a,g14.6)' ) '  Min diag/off-diag   ',ratmin
  end if

  return
end
subroutine rswrit(b1jbl,b2jbl,b3jbl,cost,e,gamt,icrys,jcrys, &
  l0,m0,nbot,p,pc,psi,rueta,ruksi,t,te,tk,tnow,u,v,vort,w,x, &
  xbot,xc,y,ybot,yc)

!*****************************************************************************80
!
!! RSWRIT writes out information which can be used for restarts or plots.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real COST, the value of the cost functional
!    associated with this solution.
!
!    Input, real E(NI,NJ), the magnetic stream function.
!
!    Input, real GAMT(NI,NJ), the diffusion coefficient.
!
!    Input, integer ( kind = 4 ) ICRYS, JCRYS, the I and J coordinates of
!    the lower left corner of the crystal.
!
!    Input, integer ( kind = 4 ) L0, the number of rows of data.
!
!    Input, integer ( kind = 4 ) M0, the number of columns of data.
!
!    Input, integer ( kind = 4 ) NBOT, the number of crucible nodes.
!
!    Input, real P(NI,NJ), the pressure.
!
!    Input, real PC(NI,NJ), the corrected pressure.
!
!    Input, real PSI(NI,NJ), the stream function.
!
!    Input, real RUETA(NI,NJ), RUKSI(NI,NJ), the
!    the momentum in the ETA and KSI directions.
!
!    Input, real T(NI,NJ), the temperature.
!
!    Input, real TE(NI,NJ), the turbulent epsilon.
!
!    Input, real TK(NI,NJ), the turbulent K.
!
!    Output, real TNOW, the current time.
!
!    Input, real U(NI,NJ), the horizontal velocity.
!
!    Input, real V(NI,NJ), the vertical velocity.
!
!    Input, real VMAG(NI,NJ), the velocity magnitude.
!
!    Input, real VORT(NI,NJ), the vorticity.
!
!    Input, real W(NI,NJ), the axial velocity.
!
!    Input, real X(NI,NJ), the X coordinates of the primary
!    nodes.
!
!    Input, real XBOT(NBOT), the X coordinates of the crucible
!    bottom nodes.
!
!    Input, real XC(NI,NJ), the X coordinates of the corner
!    nodes.
!
!    Input, real Y(NI,NJ), the Y coordinates of the primary
!    nodes.
!
!    Input, real YBOT(NBOT), the Y coordinates of the crucible
!    bottom nodes.
!
!    Input, real YC(NI,NJ), the Y coordinates of the corner
!    nodes.
!
  implicit none

  integer ( kind = 4 ) nbot
  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) b3jbl(ni)
  real ( kind = 8 ) cost
  logical, parameter :: debug = .false.
  real ( kind = 8 ) e(ni,nj)
  real ( kind = 8 ) gamt(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) m0
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) pc(ni,nj)
  real ( kind = 8 ) psi(ni,nj)
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) t(ni,nj)
  real ( kind = 8 ) te(ni,nj)
  real ( kind = 8 ) tk(ni,nj)
  real ( kind = 8 ) tnow
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) vort(ni,nj)
  real ( kind = 8 ) w(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xbot(nbot)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ybot(nbot)
  real ( kind = 8 ) yc(ni,nj)

  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) jmin
  real ( kind = 8 ) vmax
  real ( kind = 8 ) vmin

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RSWRIT:'
  write ( *, '(a)' ) '  Writing restart information.'
  write ( *, '(a)' ) ' '
!
!  Temporary check that VORT contains legitimate data.
!
  if ( debug ) then

    vmax = vort(1,1)
    imax = 1
    jmax = 1

    vmin = vort(1,1)
    imin = 1
    jmin = 1

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I, J, VORT(I,J)'
    write ( *, '(a)' ) ' '

    do i = 1,l0
      do j = 1,m0

        if ( vort(i,j) > vmax ) then
          vmax = vort(i,j)
          imax = i
          jmax = j
        end if

        if ( vort(i,j) < vmin ) then
          vmin = vort(i,j)
          imin = i
          jmin = j
        end if

      write(*,'(2i6,g14.5)')i,j,vort(i,j)

      end do
    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'RSWRIT - Temporary note:'
    write ( *, '(a,g14.6)' ) '  Minimum vorticity is ',vmin
    write ( *, '(a,i6,a,i6)' ) '  IMIN = ',imin,' JMIN=',jmin
    write ( *, '(a,g14.6)' ) '  Maximum vorticity is ',vmax
    write ( *, '(a,i6,a,i6)' ) '  IMAX = ',imax,' JMAX=',jmax
    write ( *, '(a,g14.6)' ) '  VMIN = ', vmin
    write ( *, '(a,g14.6)' ) '  VMAX = ', vmax

  end if
!
!  Delete any old copy of the restart file.
!
  open(unit = 10,file='rswrit.txt',status='old',err=10)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RSWRIT - Note:'
  write ( *, '(a)' ) '  Deleting a previous copy of the restart file.'
  close(unit = 10,status='delete')

10    continue

  open(unit = 10,file='rswrit.txt',status='unknown')

  write(10,*)cost
  write(10,*)l0
  write(10,*)jcrys
  write(10,*)icrys
  write(10,*)m0
  write(10,*)nbot
  write(10,*)tnow
!
!  NEW NEW NEW
!
  write(10,'(6g14.5)')(b1jbl(i),i = 1,nj)
  write(10,'(6g14.5)')(b2jbl(i),i = 1,nj)
  write(10,'(6g14.5)')(b3jbl(i),i = 1,ni)
  write(10,'(6g14.5)')((e(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((gamt(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((p(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((pc(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((psi(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((rueta(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((ruksi(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((t(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((te(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((tk(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((u(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((v(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((vort(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((w(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((x(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')(xbot(i),i = 1,nbot)
  write(10,'(6g14.5)')((xc(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')((y(i,j),i = 1,l0),j=1,m0)
  write(10,'(6g14.5)')(ybot(i),i = 1,nbot)
  write(10,'(6g14.5)')((yc(i,j),i = 1,l0),j=1,m0)

  close(unit = 10)

  return
end
subroutine setcst ( area, cinc, gam, icost, icrys, jcrys, l0, m0, ni, &
  nj, t, tf, u, v, vave, vort, xc, yc )

!*****************************************************************************80
!
!! SETCST computes a portion on the cost functional.
!
!  Discussion:
!
!    The cost functional is currently the integral over time and
!    space of the norm of the fluid velocity.
!
!    Currently, a crude approximation is made to the integrals.
!    For instance, the area of the flow region, LIQUID(T), is computed
!    by summing up the estimated areas of each control volume.  These
!    areas, in turn, are estimated by computing the jacobian at the
!    center of the control volume.
!
!      +-----V-----+
!      |           |
!      |           |
!      U     N     U
!      |           |
!      |           |
!      +-----V-----+
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real ( kind = 8 ) AREA(NI,NJ).
!    AREA(I,J) is the area of the control volume (I,J),
!    centered on primary node (I,J), with corner nodes
!    (I,J), (I+1,J), (I,J+1) and (I+1,J+1).
!
!    Output, real ( kind = 8 ) CINC.
!    CINC is the value of F(T) at the current time, that is, the
!    integral of the cost function over the current region.
!
!    Input, real ( kind = 8 ) GAM(NI,NJ).
!    GAM(I,J) contains the transport coefficient for the current equation
!    at primary node (I,J).  GAM(I,J) is general the fluid viscosity
!    (equations for variables 1, 2, 3, or 4), or the thermal diffusion
!    coefficient (equations for variable 5).
!    GAM is only needed for cost function ICOST = 5.  If SETCST is
!    called with ICOST = 5, then routine GAMSOR should have just
!    been called, to set the value of GAM for the heat equation.
!
!    Input, integer ( kind = 4 ) ICOST.
!    ICOST chooses the cost functional to be minimized.  The cost
!    functional is an integral over time (from TINIT to TEND) and space
!    (the region Omega).  The integration in space is normally done
!    first, yielding a function F(T) to be integrated over time.
!
!    Depending on the value of ICOST, F(T) is:
!
!      ICOST = 1: Velocity magnitude:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T))
!          (U(T,X,Y)**2 + V(T,X,Y)**2) dX dY )
!
!      ICOST = 2: Temperature deviation from TF:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T))
!          (T(T,X,Y)-TF)**2 dX dY )
!
!      ICOST = 3: Velocity magnitude deviation from Vave:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T))
!          ( Vave - SQRT(U(T,X,Y)**2 + V(T,X,Y)**2)) dX dY )
!
!      ICOST = 4: Flow velocity under the crystal:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T) and below CRYSTAL(T) )
!          U(T,X,Y)**2 + V(T,X,Y)**2 dX dY) )
!
!        For now, we take "below the crystal" to mean ALL control volumes
!        under the crystal, for I = 1 to ICRYS-1.
!
!      ICOST = 5: Heat flux from liquid into crystal and free surface:
!
!        F(T)  =  Integral ( (X,Y) in LIQUID(T) and (CRYSTAL(T) or VOID(T))
!                 GAM(X,Y) * dT(X,Y)/dNormal dX dY)
!
!      ICOST = 6: Vorticity of the flow field:
!
!        F(T)  =  SQRT( Integral ( (X,Y) in LIQUID(T))
!          VORT(T,X,Y)**2 dX dY )
!
!    Input, integer ( kind = 4 ) ICRYS.
!    ICRYS specifies the end of the crystal in the I array direction,
!    and in the vertical coordinate direction.
!
!    Input, integer ( kind = 4 ) JCRYS.
!    JCRYS specifies the end of the crystal in the J array direction,
!    and in the horizontal coordinate direction.
!
!    Input, integer ( kind = 4 ) L0.
!    L0 is the extent of the grid in the vertical or "I" coordinate.
!
!    Input, integer ( kind = 4 ) M0.
!    M0 is the extent of the grid in the horizontal or "J" coordinate.
!
!    Input, integer ( kind = 4 ) NI.
!    The maximum number of grid points in the "I" or vertical direction.
!
!    Input, integer ( kind = 4 ) NJ.
!    The maximum number of grid points in the "J" or horizontal direction.
!
!    Input, real ( kind = 8 ) T(NI,NJ).
!    The current temperature at each primary node (I,J).
!
!    Input, real ( kind = 8 ) TF.
!    The fusion or melting temperature of the crystal.
!
!    Input, real ( kind = 8 ) U(NI,NJ).
!    The U component of velocity at each primary node (I,J).
!
!    Input, real ( kind = 8 ) V(NI,NJ).
!    The V component of velocity at each primary node (I,J).
!
!    Input, real ( kind = 8 ) VAVE.
!    For cost function ICOST = 3, the desired average velocity magnitude.
!
!    Input, real ( kind = 8 ) VORT(NI,NJ).
!    The fluid flow vorticity at each primary node (I,J).
!
!    Input, real ( kind = 8 ) XC(NI,NJ).
!    XC contains the X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Input, real ( kind = 8 ) YC(NI,NJ).
!    YC contains the Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj

  real ( kind = 8 ) area(ni,nj)
  real ( kind = 8 ) cinc
  real ( kind = 8 ) dist
  real ( kind = 8 ) dtdn
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) gamma
  real ( kind = 8 ) hf
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icost
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) m0
  real ( kind = 8 ) t(ni,nj)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) tf
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) vave
  real ( kind = 8 ) vort(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
!
!  #1: Velocity magnitude:
!
  if ( icost == 1 ) then

    cinc = 0.0D+00
    do i = 1,icrys-1
      do j = 1,m0
        cinc = cinc+(u(i,j)**2+v(i,j)**2)*area(i,j)
      end do
    end do

    cinc = sqrt(cinc)
!
!  #2: Temperature deviation from TF:
!
  else if ( icost == 2 ) then

    cinc = 0.0D+00
    do i = 1,icrys-1
      do j = 1,m0
        cinc = cinc+(t(i,j)-tf)**2*area(i,j)
      end do
    end do

    cinc = sqrt(cinc)
!
!  #3: Velocity magnitude deviation from VAVE:
!
  else if ( icost == 3 ) then

    cinc = 0.0D+00
    do i = 1,icrys-1
      do j = 1,m0
        cinc = cinc+area(i,j)*(vave-sqrt(u(i,j)**2+v(i,j)**2))**2
      end do
    end do

    cinc = sqrt(cinc)
!
!  #4: Flow velocity under the crystal:
!
!  For now, we take "under the crystal" to mean ALL control volumes
!  under the crystal, for I = 1 to ICRYS-1.
!
  else if ( icost == 4 ) then

    cinc = 0.0D+00
    do i = 1,icrys-1
      do j = 1,jcrys
        cinc = cinc+area(i,j)*((u(i,j)**2+v(i,j)**2))
      end do
    end do

    cinc = sqrt(cinc)
!
!  #5: Heat flux through crystal and free surface:
!
  else if ( icost == 5 ) then

    cinc = 0.0D+00
!
!  Add up the heat flux passing into the BOTTOM of the crystal.
!  I believe this involves heat passing from Control Volume
!  CV(ICRYS-1,J) to Control Volume CV(ICRYS,J), for J = 2 to JCRYS.
!
    do j = 2,jcrys

      t1 = 0.5D+00
      t2 = 0.5D+00
      gamma = (t1*gam(icrys,j+1)+t2*gam(icrys,j))/(t1+t2)

      dist = sqrt((xc(icrys,j+1)-xc(icrys,j))**2+(yc(icrys,j+1)-yc(icrys,j))**2)

      dtdn = (t(icrys,j+1)-t(icrys,j))/dist

      hf = gamma*dtdn

      cinc = cinc+hf*(xc(icrys,j+1)-xc(icrys,j))

    end do
!
!  Add up the heat flux passing into the void from the free surface.
!  This involves heat passing from Control Volume CV(ICRYS-1,J)
!  to Control Volume CV(ICRYS,J) for J = JCRYS+1 to M0.
!
    do j = jcrys+1,m0

      t1 = 0.5D+00
      t2 = 0.5D+00
      gamma = (t1*gam(icrys,j+1)+t2*gam(icrys,j))/(t1+t2)

      dist = sqrt((xc(icrys,j+1)-xc(icrys,j))**2+(yc(icrys,j+1)-yc(icrys,j))**2)

      dtdn = (t(icrys,j+1)-t(icrys,j))/dist

      hf = gamma*dtdn

      cinc = cinc+hf*(xc(icrys,j+1)-xc(icrys,j))

    end do
!
!  #6: Vorticity of the flow field:
!
  else if ( icost == 6 ) then

    cinc = 0.0D+00
    do i = 1,l0
      do j = 1,m0
        cinc = cinc+area(i,j)*vort(i,j)**2
      end do
    end do

    cinc = sqrt(cinc)

  end if

  return
end
subroutine setgeo ( ae1, ae2, ak1, ak2, heta, hksi, l1, m1, mode, ni, nj, &
  r, vol, x, xc, y, yc )

!*****************************************************************************80
!
!! SETGEO calculates various geometric quantities.
!
!  Discussion:
!
!    SETGEO is given (XC,YC), the locations of the "corners" of the
!    control volumes, and calculates various related geometric parameters,
!    including:
!
!      X, Y the position of the primary nodes,
!      HETA and HKSI, dH/dETA and dH/dKSI,
!      the Jacobian,
!      AK1, AE1, dALPHA/dKSI, dALPHA/dETA,
!      AK2, AE2, dBETA/dKSI, dBETA/dETA,
!      R, a weight factor for cylindrical geometries.
!      VOL, the volume of the control volumes.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Output, real ( kind = 8 ) AE1(NI,NJ).
!    d ALPHA/d ETA in thesis.
!    I think AE1 is the value of d ETA/d X or d ETA/d Y.
!
!    Output, real ( kind = 8 ) AE2(NI,NJ).
!    d BETA/d ETA in thesis.
!    I think AE2 is the value of d KSI/d X or d KSI/d Y.
!
!    Output, real ( kind = 8 ) AK1(NI,NJ).
!    d ALPHA/d KSI in thesis.
!    I think AK1 is the value of d ETA/d X or d ETA/d Y.
!
!    Output, real ( kind = 8 ) AK2(NI,NJ).
!    d BETA/d KSI in thesis.
!    I think AK2 is the value of d KSI/d X or d KSI/d Y.
!
!    Output, real ( kind = 8 ) HETA(NI,NJ).
!    HETA(I,J) is the physical (that is, "X,Y") length of the control
!    volume in the ETA direction.
!
!    HETA(I,J) is the length of the line connecting the mid-side
!    nodes that separate primary node (I,J) from primary nodes
!    (I,J-1) and (I,J+1).
!
!    Warning: The code only defines HETA for the solid or liquid
!    region it is working on, and sets it to zero elsewhere!
!
!  HKSI   Output, real ( kind = 8 ) HKSI(NI,NJ).
!         HKSI(I,J) is the physical (that is, "X,Y") length of the control
!         volume in the KSI direction.
!
!         HKSI(I,J) is the length of the line connecting the mid-side
!         nodes that separate primary node (I,J) from primary nodes
!         (I-1,J) and (I+1,J).
!
!         Warning: The code only defines HKSI for the solid or liquid
!         region it is working on, and sets it to zero elsewhere!
!
!  L1     Input, integer ( kind = 4 ) L1.
!         L1 is some portion of the vertical grid size.
!         L1 could equal L0, or ICRYS, depending on which subproblem is
!         being solved.
!
!  M1     Input, integer ( kind = 4 ) M1.
!         M1 is some portion of the horizontal grid size.
!         M1 could equal M0, or JCRYS, depending on which subproblem is
!         being solved.
!
!  MODE   Input, integer ( kind = 4 ) MODE.
!         MODE determines the kind of 2D geometry to be used.
!         If MODE = 0, rectangular geometry is used;
!         If MODE = 1, axisymmetric cylindrical geometry is used.
!
!  NI     Input, integer ( kind = 4 ) NI.
!         NI is the maximum number of grid points in the I or row or vertical
!         direction.
!
!  NJ     Input, integer ( kind = 4 ) NJ.
!         NJ is the maximum number of grid points in the J or column or
!         horizontal direction.
!
!  R      Output, real ( kind = 8 ) R(NI,NJ).
!         R is used in axisymmetric problems to store the radial distance.
!         For MODE = 0, R(I,J) = 1.
!         For MODE = 1, R(I,J) = Y(I,J), the radial distance.
!
!  VOL    Output, real ( kind = 8 ) VOL(NI,NJ).
!         VOL contains the volume of each control volume.
!         Unfortunately, this is not completely accurate.
!
!         The code only defines the volume for the control volumes that happen
!         to be in the current subregion of interest.  Outside of that
!         region, the value of VOL cannot be relied on.  The subregion of
!         interest is specified by the values of L1 and M1.
!
!         Also, for an axisymmetric problem (MODE = 1), VOL is automatically
!         multiplied by a factor representing the radial distance.
!
!  X      Output, real ( kind = 8 ) X(NI,NJ).
!         X contains the X coordinates of the primary nodes, which are
!         the centers of the control volumes.
!
!         In most cases, X(I,J) is the average of the four enclosing
!         corner values, XC(I,J), XC(I+1,J), XC(I,J+1), and XC(I+1,J+1).
!
!  XC     Input, real ( kind = 8 ) XC(NI,NJ).
!         XC contains the X coordinate of nodes which are the "corners" of
!         control volumes.
!
!  Y      Output, real ( kind = 8 ) Y(NI,NJ).
!         Y contains the Y coordinates of primary nodes, which are
!         the centers of the control volumes.
!
!         In most cases, Y(I,J) is the average of the four enclosing
!         corner values, YC(I,J), YC(I+1,J), YC(I,J+1), and YC(I+1,J+1).
!
!  YC     Input, real ( kind = 8 ) YC(NI,NJ).
!         YC contains the Y coordinate of nodes which are the "corners" of
!         control volumes.
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj

  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdksi
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydksi
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mode
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xd
  real ( kind = 8 ) xjacb
  real ( kind = 8 ) xm
  real ( kind = 8 ) xp
  real ( kind = 8 ) xu
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yd
  real ( kind = 8 ) ym
  real ( kind = 8 ) yp
  real ( kind = 8 ) yu
!
!  Calculate HKSI and HETA, essentially the width and height of
!  the control volume, and VOL, the volume of the control volume.
!
  do i = 2,l1-1
    do j = 2,m1-1
!
!  Compute the coordinates of the mid-side nodes.
!
!    ^
!    |     [I,J+1]------Up--------[I+1,J+1]
!    |       |                        |
!    E       |                        |
!    T     Minus        (I,J)       Plus
!    A       |                        |
!    |       |                        |
!    |     [I,J]--------Down-------[I+1,J]
!    |
!    |
!    +-------------XSI---->
!
      xp = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j))
      yp = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j))

      xm = 0.5D+00*(xc(i,j+1)+xc(i,j))
      ym = 0.5D+00*(yc(i,j+1)+yc(i,j))

      xu = 0.5D+00*(xc(i,j+1)+xc(i+1,j+1))
      yu = 0.5D+00*(yc(i,j+1)+yc(i+1,j+1))

      xd = 0.5D+00*(xc(i,j)+xc(i+1,j))
      yd = 0.5D+00*(yc(i,j)+yc(i+1,j))
!
!  The mid-side nodes differ by Delta KSI or Delta ETA  =  1.
!  Use finite differences to estimate derivatives like dX/dKSI.
!
      dxdksi = xp-xm
      dxdeta = xu-xd
      dydksi = yp-ym
      dydeta = yu-yd
!
!  Now compute the length of the line segments that connect
!  opposing mid-side nodes.
!
      hksi(i,j) = sqrt(dxdksi**2+dydksi**2)
      heta(i,j) = sqrt(dxdeta**2+dydeta**2)
!
!  Now compute the area of the control volume, which is just
!  the determinant of the Jacobian matrix.
!
      vol(i,j) = dxdksi*dydeta-dxdeta*dydksi

    end do
  end do
!
!  Take care of values along the borders J = 1 and J=M1.
!
  do i = 1,l1

    heta(i,1) = 0.0D+00
    heta(i,m1) = 0.0D+00
    vol(i,1) = 0.0D+00
    vol(i,m1) = 0.0D+00

    if ( i /= 1 .and. i.ne.l1 ) then
      hksi(i,1) = sqrt((xc(i+1,2)-xc(i,2))**2+(yc(i+1,2)-yc(i,2))**2)
      hksi(i,m1) = sqrt((xc(i+1,m1)-xc(i,m1))**2+(yc(i+1,m1)-yc(i,m1))**2)
    else
      hksi(i,1) = 0.0D+00
      hksi(i,m1) = 0.0D+00
    end if

  end do
!
!  Take care of values along the borders I = 1 and I=L1.
!
  do j = 1,m1

    hksi(1,j) = 0.0D+00
    hksi(l1,j) = 0.0D+00
    vol(1,j) = 0.0D+00
    vol(l1,j) = 0.0D+00

    if ( j /= 1 .and. j.ne.m1 ) then
      heta(1,j) = sqrt((xc(2,j)-xc(2,j+1))**2+(yc(2,j)-yc(2,j+1))**2)
      heta(l1,j) = sqrt((xc(l1,j)-xc(l1,j+1))**2+(yc(l1,j)-yc(l1,j+1))**2)
    else
      heta(1,j) = 0.0D+00
      heta(l1,j) = 0.0D+00
    end if

  end do
!
!  Calculate AK1 and AK2, whose meaning has not been revealed to me.
!
!
!    ^
!    |                      C[I,J+1]
!    |                         |
!    E                         |
!    T       P(I-1,J)          $     P(I,J)
!    A                         |
!    |                         |
!    |                      C[I,J]
!    |
!    |
!    +---------XSI--------->
!
  do j = 2,m1-1
    do i = 2,l1

      dxdeta = xc(i,j+1)-xc(i,j)
      dydeta = yc(i,j+1)-yc(i,j)
      dxdksi = x(i,j)-x(i-1,j)
      dydksi = y(i,j)-y(i-1,j)

      if ( i == 2 .or. i == l1 ) then
        dxdksi = dxdksi*2.0D+00
        dydksi = dydksi*2.0D+00
      end if

      t1 = sqrt(dxdeta**2+dydeta**2)
      t2 = sqrt(dxdksi**2+dydksi**2)
      t3 = dxdksi*dxdeta+dydeta*dydksi

      xjacb = dxdksi*dydeta-dxdeta*dydksi

      if ( xjacb /= 0.0D+00 ) then
        ak1(i,j) = t2*t1*t1/xjacb
        ak2(i,j) = t1*t3/xjacb
      else
        ak1(i,j) = 0.0D+00
        ak2(i,j) = 0.0D+00
      end if

    end do
  end do
!
!  Calculate AE1 and AE2, whose meaning has not been revealed to me.
!
!    ^
!    |           P(I,J)
!    |
!    E
!    T   C[I,J]-----$-------C[I+1,J]
!    A
!    |
!    |           P(I,J-1)
!    |
!    +-------------XSI------------->
!
  do i = 2,l1-1
    do j = 2,m1

      dxdeta = x(i,j)-x(i,j-1)
      dydeta = y(i,j)-y(i,j-1)
      dxdksi = xc(i+1,j)-xc(i,j)
      dydksi = yc(i+1,j)-yc(i,j)

      if ( j == 2 .or. j == m1 ) then
        dxdeta = dxdeta*2.0D+00
        dydeta = dydeta*2.0D+00
      end if

      t1 = sqrt(dxdeta**2+dydeta**2)
      t2 = sqrt(dxdksi**2+dydksi**2)
      t3 = dxdksi*dxdeta+dydeta*dydksi

      xjacb = dxdksi*dydeta-dxdeta*dydksi

      if ( xjacb /= 0.0D+00 ) then
        ae1(i,j) = t1*t2*t2/xjacb
        ae2(i,j) = t2*t3/xjacb
      else
        ae1(i,j) = 0.0D+00
        ae2(i,j) = 0.0D+00
      end if

    end do
  end do
!
!  Set the cylindrical geometry weight factor.
!
  if ( mode == 0 ) then

    do i = 1,l1
      do j = 1,m1
        r(i,j) = 1.0D+00
      end do
    end do

  else if ( mode == 1 ) then

    do i = 1,l1
      do j = 1,m1

        r(i,j) = y(i,j)

        ak1(i,j) = ak1(i,j)*r(i,j)
        ak2(i,j) = ak2(i,j)*r(i,j)
        ae1(i,j) = ae1(i,j)*r(i,j)
        ae2(i,j) = ae2(i,j)*r(i,j)
        vol(i,j) = vol(i,j)*r(i,j)

      end do
    end do

  end if

  return
end
subroutine setp ( fjeta, fjksi, heta, hksi, l1, m1, ni, nj, p )

!*****************************************************************************80
!
!! SETP estimates the values of pressure at the interfaces.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj

  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) tt

  do i = 2,l1
    do j = 1,m1

      if ( hksi(i,j)+hksi(i-1,j) /= 0.0D+00 ) then
        tt = hksi(i,j)/(hksi(i,j)+hksi(i-1,j))
      else
        tt = 0.0D+00
      end if

      fjksi(i,j) = tt*p(i-1,j)+(1.0D+00-tt)*p(i,j)

    end do
  end do

  do i = 1,l1
    do j = 2,m1

      if ( heta(i,j)+heta(i,j-1) /= 0.0D+00 ) then
        tt = heta(i,j)/(heta(i,j)+heta(i,j-1))
      else
        tt = 0.0D+00
      end if

      fjeta(i,j) = tt*p(i,j-1)+(1.0D+00-tt)*p(i,j)

    end do
  end do

  return
end
subroutine setru ( heta, hksi, l1, m1, ni, nj, rho, rueta, ruksi, u, &
  v, x, y )

!*****************************************************************************80
!
!! SETRU estimates the normal mass velocity component at interfaces.
!
!  Discussion:
!
!    SETRU estimates the normal mass velocity component at the east-west
!    and north-south interfaces of the control volumes by interpolating
!    values of velocity and density computed for the primary nodes.
!
!    These values are needed only for cells which lie within the melt
!    region.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Double precision HETA(NI,NJ).
!    HETA(I,J) is the physical (that is, "X,Y") length of the control
!    volume in the ETA direction.
!
!    HETA(I,J) is the length of the line connecting the mid-side
!    nodes that separate primary node (I,J) from primary nodes
!    (I,J-1) and (I,J+1).
!
!    Warning: The code only defines HETA for the solid or liquid
!    region it is working on, and sets it to zero elsewhere!
!
!    Double precision HKSI(NI,NJ).
!    HKSI(I,J) is the physical (that is, "X,Y") length of the control
!    volume in the XI direction.
!
!    HKSI(I,J) is the length of the line connecting the mid-side
!    nodes that separate primary node (I,J) from primary nodes
!    (I-1,J) and (I+1,J).
!    Warning: The code only defines HKSI for the solid or liquid
!    region it is working on, and sets it to zero elsewhere!
!
!    Integer L1.
!    L1 is some portion of the vertical grid size.
!    L1 could equal L0, or ICRYS, depending on which subproblem is
!    being solved.
!
!    Integer M1.
!    M1 is some portion of the horizontal grid size.
!    M1 could equal M0, or JCRYS, depending on which subproblem is
!    being solved.
!
!    Integer NI.
!    NI is the maximum number of grid points in the I or row or vertical
!    or XI direction.
!
!    Integer NJ.
!    NJ is the maximum number of grid points in the J or column or
!    horizontal or ETA direction.
!
!    Double precision RHO(NI,NJ).
!    RHO contains the fluid density at each primary node.
!    In version 1.0 of MASTRAPP, RHO is has a constant uniform
!    value, namely RHOCON.
!
!    Double precision RUETA(NI,NJ).
!    RUETA is the component of the mass velocity interpolated
!    at the interface between primary nodes (I,J-1) and (I,J),
!    and normal to the cell interface.
!
!    Double precision RUKSI(NI,NJ).
!    RUKSI is the component of the mass velocity interpolated
!    at the interface between primary nodes (I-1,J) and (I,J),
!    and normal to the cell interface.
!
!    Double precision U(NI,NJ).
!    U contains the horizontal (in the XY coordinate system) component of
!    velocity at each primary node (I,J).
!
!    Double precision V(NI,NJ).
!    V contains the vertical (in the XY coordinate system) component of
!    velocity at each primary node (I,J).
!
!    Double precision X(NI,NJ).
!    X contains the X coordinates of the primary nodes, which are
!    the centers of the control volumes.
!    In most cases, X(I,J) is the average of the four enclosing
!    corner values, XC(I,J), XC(I+1,J), XC(I,J+1), and XC(I+1,J+1).
!
!    Double precision Y(NI,NJ).
!    Y contains the Y coordinates of primary nodes, which are
!    the centers of the control volumes.
!    In most cases, Y(I,J) is the average of the four enclosing
!    corner values, YC(I,J), YC(I+1,J), YC(I,J+1), and YC(I+1,J+1).
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj

  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhov
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) tt
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) ucs
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) vcs
  real ( kind = 8 ) vlen
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xcomp
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) ycomp
!
!  Compute RUKSI, the normal velocity at east-west cell interfaces.
!
  do i = 2,l1
    do j = 1,m1
!
!  To interpolate a value at the cell interface, we interpolate various
!  quantities using the relative distance of the (I-1,J) and (I,J)
!  cells from their common interface.
!
      if ( hksi(i,j)+hksi(i-1,j) /= 0.0D+00 ) then
        tt = hksi(i,j)/(hksi(i,j)+hksi(i-1,j))
      else
        tt = 0.0D+00
      end if
!
!  Interpolate the values of the velocity and density from the
!  centers of the two cells to the intermediate cell interface.
!
      ucs =  tt*u(i-1,j)  +(1.0D+00-tt)*u(i,j)
      vcs =  tt*v(i-1,j)  +(1.0D+00-tt)*v(i,j)
      rhov = tt*rho(i-1,j)+(1.0D+00-tt)*rho(i,j)
!
!  Compute the unit vector which is normal to the cell interface.
!
      xcomp = x(i,j)-x(i-1,j)
      ycomp = y(i,j)-y(i-1,j)
      vlen = sqrt(xcomp**2+ycomp**2)
      if ( vlen > 0.0D+00 ) then
        xcomp = xcomp/vlen
        ycomp = ycomp/vlen
      end if
!
!  Interpolate the normal mass velocity to the cell interface.
!
      ruksi(i,j) = rhov*(ucs*xcomp+vcs*ycomp)

    end do
  end do
!
!  Compute RUETA, the normal velocity at north-south cell interfaces.
!
  do i = 1,l1
    do j = 2,m1
!
!  To interpolate a value at the cell interface, we interpolate various
!  quantities using the relative distance of the (I,J-1) and (I,J)
!  cells from their common interface.
!
      if ( heta(i,j)+heta(i,j-1) /= 0.0D+00 ) then
        tt = heta(i,j)/(heta(i,j)+heta(i,j-1))
      else
        tt = 0.0D+00
      end if
!
!  Interpolate the values of the velocity and density from the
!  centers of the two cells to the intermediate cell interface.
!
      ucs = tt*u(i,j-1)+(1.0D+00-tt)*u(i,j)
      vcs = tt*v(i,j-1)+(1.0D+00-tt)*v(i,j)
      rhov = tt*rho(i,j-1)+(1.0D+00-tt)*rho(i,j)
!
!  Compute the unit vector which is normal to the cell interface.
!
      xcomp = x(i,j)-x(i,j-1)
      ycomp = y(i,j)-y(i,j-1)
      vlen = sqrt(xcomp**2+ycomp**2)
      if ( vlen > 0.0D+00 ) then
        xcomp = xcomp/vlen
        ycomp = ycomp/vlen
      end if
!
!  Interpolate the normal mass velocity to the cell interface.
!
      rueta(i,j) = rhov*(ucs*xcomp+vcs*ycomp)

    end do
  end do

  return
end
subroutine setup(ae1,ae2,ak1,ak2,b1jbl,b2jbl,b3jbl,birad, &
  cappa, &
  cd,ce1,ce2,cmu,epsil,ewall,f,fcsl,fjeta,fjksi,fksl,fma, &
  fmax,fn,fo,frsl,gamt,grash,hamag,heta,hksi,inturb,icrys, &
  iter,izone,jcrys,l1,lblk,lconv,lortho,lsolve,m1, &
  mode,nf,np,npc,nsolve,ntimes,pr,r,rdtm,re,recb,rect,res, &
  relax,rho,rhocon,rpr,rueta,ruksi,sige,sigk,sigt,smax,ssum, &
  stel,stes,tal,tas,tf,tw,vol,x,xc,y,yc)

!*****************************************************************************80
!
!! SETUP calculates the coefficients of the equations, and solves them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real ( kind = 8 ) XC(NI,NJ).
!    The X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Input, real ( kind = 8 ) YC(NI,NJ).
!    The Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: nk = 14
  integer ( kind = 4 ), parameter :: nmaxij = 64
  integer ( kind = 4 ), parameter :: ns = 10

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) acof
  real ( kind = 8 ) ae1(ni,nj)
  real ( kind = 8 ) ae2(ni,nj)
  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ak1(ni,nj)
  real ( kind = 8 ) ak2(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) apr
  real ( kind = 8 ) b1jbl(nj)
  real ( kind = 8 ) b2jbl(nj)
  real ( kind = 8 ) b3jbl(ni)
  real ( kind = 8 ) birad
  real ( kind = 8 ) biu(nj)
  real ( kind = 8 ) biv(nj)
  real ( kind = 8 ) bju(ni)
  real ( kind = 8 ) bjv(ni)
  real ( kind = 8 ) blu(nj)
  real ( kind = 8 ) blv(nj)
  real ( kind = 8 ) bmu(ni)
  real ( kind = 8 ) bmv(ni)
  real ( kind = 8 ) bpi
  real ( kind = 8 ) bpi1
  real ( kind = 8 ) bpj
  real ( kind = 8 ) bpj1
  real ( kind = 8 ) bpm
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) ce1
  real ( kind = 8 ) ce2
  real ( kind = 8 ) cmu
  real ( kind = 8 ) cofu(ni,nj,5)
  real ( kind = 8 ) cofv(ni,nj,5)
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) conu(ni,nj)
  real ( kind = 8 ) conv(ni,nj)
  real ( kind = 8 ) denom
  real ( kind = 8 ) diff
  real ( kind = 8 ) dpeta
  real ( kind = 8 ) dpksi
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdksi
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydksi
  real ( kind = 8 ) em
  real ( kind = 8 ) ep
  real ( kind = 8 ) epsil
  real ( kind = 8 ) ewall
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fcsl
  real ( kind = 8 ) fjeta(ni,nj)
  real ( kind = 8 ) fjksi(ni,nj)
  real ( kind = 8 ) fksl
  real ( kind = 8 ) flow
  real ( kind = 8 ) fma
  real ( kind = 8 ) fmax(ns)
  real ( kind = 8 ) fn(ni,nj,ns)
  real ( kind = 8 ) fo(ni,ni,ns)
  real ( kind = 8 ) frc
  real ( kind = 8 ) frsl
  real ( kind = 8 ) gam(ni,nj)
  real ( kind = 8 ) gamt(ni,nj)
  real ( kind = 8 ) grash
  real ( kind = 8 ) hamag
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hetap(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  real ( kind = 8 ) hksip(ni,nj)
  real ( kind = 8 ) hutop(ni,nj)
  real ( kind = 8 ) hvtop(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) inturb
  integer ( kind = 4 ) isol
  integer ( kind = 4 ) iter
  integer ( kind = 4 ) izone
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l1
  logical lblk(ns)
  logical lconv
  logical lortho
  logical lsolve(ns)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) mode
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npc
  integer ( kind = 4 ) nsolve(ns)
  integer ( kind = 4 ) ntimes(ns)
  real ( kind = 8 ) peta(ni,nj)
  real ( kind = 8 ) pksi(ni,nj)
  real ( kind = 8 ) pr
  real ( kind = 8 ) qt(nmaxij)
  real ( kind = 8 ) r(ni,nj)
  real ( kind = 8 ) rdtm
  real ( kind = 8 ) re
  real ( kind = 8 ) recb
  real ( kind = 8 ) rect
  real ( kind = 8 ) relax(nk)
  real ( kind = 8 ) res(ns)
  real ( kind = 8 ) rho(ni,nj)
  real ( kind = 8 ) rhocon
  real ( kind = 8 ) rpr
  real ( kind = 8 ) rueij
  real ( kind = 8 ) rueij1
  real ( kind = 8 ) ruet
  real ( kind = 8 ) rueta(ni,nj)
  real ( kind = 8 ) ruki1j
  real ( kind = 8 ) rukij
  real ( kind = 8 ) ruksi(ni,nj)
  real ( kind = 8 ) rukt
  real ( kind = 8 ) sige
  real ( kind = 8 ) sigk
  real ( kind = 8 ) sigt
  real ( kind = 8 ) smax
  real ( kind = 8 ) soor
  real ( kind = 8 ) ssum
  real ( kind = 8 ) stel
  real ( kind = 8 ) stes
  real ( kind = 8 ) t1
  real ( kind = 8 ) t2
  real ( kind = 8 ) t3
  real ( kind = 8 ) t4
  real ( kind = 8 ) tal
  real ( kind = 8 ) tas
  real ( kind = 8 ) tem1
  real ( kind = 8 ) tem2
  real ( kind = 8 ) temp
  real ( kind = 8 ) tf
  real ( kind = 8 ) tmp(ni,nj)
  real ( kind = 8 ) tw
  real ( kind = 8 ) vol(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) xd
  real ( kind = 8 ) xl
  real ( kind = 8 ) xr
  real ( kind = 8 ) xu
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
  real ( kind = 8 ) yd
  real ( kind = 8 ) yl
  real ( kind = 8 ) yr
  real ( kind = 8 ) yu

  do n = 1,ns

    if ( lsolve(n) ) then

      nf = n
!
!  N = 1, U VELOCITY.
!
      if ( n == 1 ) then

        call gamsor(ap,birad,ce1,ce2,cmu,con,ewall,f,fcsl,fksl, &
          fma,fo,frsl,gam,gamt,grash,hamag,heta,hksi,icrys,inturb, &
          izone,jcrys,l1,lsolve,m1,mode,nf,pr,r,rdtm,re,recb,rect, &
          rho,rhocon,rpr,sige,sigk,sigt,stel,stes,tal,tas,tf,tw,x, &
          xc,y,yc)

        do j = 2,m1-1

          diff = 2.0D+00*gam(2,j)/hksi(2,j)
          flow = -ruksi(2,j)
          call diflow(acof,diff,flow)
          aim(2,j) = acof*ak1(2,j)

          diff = 2.0D+00*gam(l1-1,j)/hksi(l1-1,j)
          flow = ruksi(l1,j)
          call diflow(acof,diff,flow)
          aip(l1-1,j) = acof*ak1(l1,j)

          do i = 2,l1-2
            diff = gam(i,j)*gam(i+1,j) / &
              (0.5D+00*hksi(i,j)*gam(i+1,j)+0.5D+00*hksi(i+1,j)*gam(i,j))
            flow = ruksi(i+1,j)
            call diflow(acof,diff,flow)
            aip(i,j) = acof*ak1(i+1,j)
            aim(i+1,j) = (acof+ruksi(i+1,j))*ak1(i+1,j)
          end do

        end do

        do i = 2,l1-1

          diff = 2.0D+00*gam(i,2)/heta(i,2)
          flow = -rueta(i,2)
          call diflow(acof,diff,flow)
          ajm(i,2) = acof*ae1(i,2)

          diff = 2.0D+00*gam(i,m1-1)/heta(i,m1-1)
          flow = rueta(i,m1)
          call diflow(acof,diff,flow)
          ajp(i,m1-1) = acof*ae1(i,m1)

          do j = 2,m1-2
            diff = gam(i,j)*gam(i,j+1)/ &
              (0.5D+00*heta(i,j)*gam(i,j+1)+0.5D+00*heta(i,j+1)*gam(i,j))
            flow = rueta(i,j+1)
            call diflow(acof,diff,flow)
            ajp(i,j) = acof*ae1(i,j+1)
            ajm(i,j+1) = (acof+rueta(i,j+1))*ae1(i,j+1)
          end do

        end do

        do j = 1, m1
          biu(j) = gam(1,j)
          blu(j) = gam(l1,j)
        end do

        do i = 1, l1
          bju(i) = gam(i,1)
          bmu(i) = gam(i,m1)
        end do

        do i = 1, l1
          do j = 1, m1
            conu(i,j) = con(i,j)
            cofu(i,j,1) = ap(i,j)
            cofu(i,j,2) = aip(i,j)
            cofu(i,j,3) = aim(i,j)
            cofu(i,j,4) = ajp(i,j)
            cofu(i,j,5) = ajm(i,j)
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,1))*vol(i,j)
          end do
        end do

        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          b3jbl,cappa,cd,cmu, &
          con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        do i = 2,l1-1
          do j = 2,m1-1
            hutop(i,j) = (aip(i,j)*f(i+1,j,1)+aim(i,j)*f(i-1,j,1) &
              +ajp(i,j)*f(i,j+1,1)+ajm(i,j)*f(i,j-1,1)+con(i,j))/vol(i,j)
          end do
        end do
!
!  N = 2, V VELOCITY.
!
      else if ( n == 2 ) then

        call gamsor(ap,birad,ce1,ce2,cmu,con,ewall,f,fcsl,fksl, &
          fma,fo,frsl,gam,gamt,grash,hamag,heta,hksi,icrys,inturb, &
          izone,jcrys,l1,lsolve,m1,mode,nf,pr,r,rdtm,re,recb,rect, &
          rho,rhocon,rpr,sige,sigk,sigt,stel,stes,tal,tas,tf,tw,x, &
          xc,y,yc)

        do j = 2,m1-1

          diff = 2.0D+00*gam(2,j)/hksi(2,j)
          flow = -ruksi(2,j)
          call diflow(acof,diff,flow)
          aim(2,j) = acof*ak1(2,j)

          diff = 2.0D+00*gam(l1-1,j)/hksi(l1-1,j)
          flow = ruksi(l1,j)
          call diflow(acof,diff,flow)
          aip(l1-1,j) = acof*ak1(l1,j)

          do i = 2,l1-2
            diff = gam(i,j)*gam(i+1,j)/(0.5D+00*hksi(i,j)*gam(i+1,j) &
              +0.5D+00*hksi(i+1,j)*gam(i,j))
            flow = ruksi(i+1,j)
            call diflow(acof,diff,flow)
            aip(i,j) = acof*ak1(i+1,j)
            aim(i+1,j) = (acof+ruksi(i+1,j))*ak1(i+1,j)
          end do

        end do

        do i = 2,l1-1

          diff = 2.0D+00*gam(i,2)/heta(i,2)
          flow = -rueta(i,2)
          call diflow(acof,diff,flow)
          ajm(i,2) = acof*ae1(i,2)

          diff = 2.0D+00*gam(i,m1-1)/heta(i,m1-1)
          flow = rueta(i,m1)
          call diflow(acof,diff,flow)
          ajp(i,m1-1) = acof*ae1(i,m1)

          do j = 2,m1-2
            diff = gam(i,j)*gam(i,j+1)/(0.5D+00*heta(i,j)*gam(i,j+1) &
              +0.5D+00*heta(i,j+1)*gam(i,j))
            flow = rueta(i,j+1)
            call diflow(acof,diff,flow)
            ajp(i,j) = acof*ae1(i,j+1)
            ajm(i,j+1) = (acof+rueta(i,j+1))*ae1(i,j+1)
          end do

        end do

        do j = 1,m1
          biv(j) = gam(1,j)
          blv(j) = gam(l1,j)
        end do

        do i = 1,l1
          bjv(i) = gam(i,1)
          bmv(i) = gam(i,m1)
        end do

        do i = 1,l1
          do j = 1,m1
            conv(i,j) = con(i,j)
            cofv(i,j,1) = ap(i,j)
            cofv(i,j,2) = aip(i,j)
            cofv(i,j,3) = aim(i,j)
            cofv(i,j,4) = ajp(i,j)
            cofv(i,j,5) = ajm(i,j)
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,2))*vol(i,j)
          end do
        end do

        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          b3jbl,cappa,cd,cmu, &
          con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        do i = 2,l1-1
          do j = 2,m1-1
            hvtop(i,j) = (aip(i,j)*f(i+1,j,2)+aim(i,j)*f(i-1,j,2) &
              +ajp(i,j)*f(i,j+1,2)+ajm(i,j)*f(i,j-1,2) &
              +con(i,j))/vol(i,j)
            dxdksi = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            hksip(i,j) = 0.5D+00*(dxdksi*hutop(i,j)+dydksi*hvtop(i,j))
            dxdeta = 0.5D+00*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5D+00*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            hetap(i,j) = 0.5D+00*(dxdeta*hutop(i,j)+dydeta*hvtop(i,j))
          end do
        end do
!
!  N = 3, P PRESSURE.
!
!  HKSIP (stored in CON0), HETAP (stored in AP0)
!  are based on momentum interpolation.
!  T3 is based on the two-dimenisonal correction of MIS
!  T4 is based on the secondary correction of 2-d of MIS
!
      else if ( n == 3 ) then

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = 0.0D+00
            ap(i,j) = (aip(i,j)+aim(i,j)+ajm(i,j)+ajp(i,j))/vol(i,j)/rho(i,j)
            tmp(i,j) = ap(i,j)
          end do
        end do

        if ( .not.lortho ) then

          do j = 2,m1-1

            do i = 1,l1
              qt(i) = 0.5D+00*(rueta(i,j)+rueta(i,j+1))
            end do

            do i = 2,l1-1
              tem1 = hksi(i+1,j)+hksi(i,j)
              tem2 = hksi(i-1,j)+hksi(i,j)
              t1 = hksi(i,j)/tem1
              t2 = hksi(i+1,j)/tem1
              t3 = hksi(i,j)/tem2
              t4 = hksi(i-1,j)/tem2
              con(i,j) = con(i,j)+(t1*qt(i+1)+t2*qt(i))*ak2(i+1,j) &
                -(t4*qt(i)+t3*qt(i-1))*ak2(i,j)
            end do

          end do

          do i = 2,l1-1

            do j = 1,m1
              qt(j) = 0.5D+00*(ruksi(i,j)+ruksi(i+1,j))
            end do

            do j = 2,m1-1
              tem1 = heta(i,j+1)+heta(i,j)
              tem2 = heta(i,j-1)+heta(i,j)
              t1 = heta(i,j)/tem1
              t2 = heta(i,j+1)/tem1
              t3 = heta(i,j)/tem2
              t4 = heta(i,j-1)/tem2
              con(i,j) = con(i,j)+(t1*qt(j+1)+t2*qt(j))*ae2(i,j+1) &
                -(t4*qt(j)+t3*qt(j-1))*ae2(i,j)
            end do
          end do

        end if

        do j = 2,m1-1

          temp = 0.0D+00

          do i = 2,l1-2

            a1 = 0.5D+00*(ak1(i,j)+ak1(i+1,j))
            a2 = 0.5D+00*(ak1(i+1,j)+ak1(i+2,j))
            denom = 0.5D+00*(ap(i,j)*hksi(i,j)/a1+ap(i+1,j)*hksi(i+1,j)/a2)
            apr = 2.0D+00*denom/(hksi(i,j)+hksi(i+1,j))

            dxdksi = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            rukij = (dxdksi*f(i,j,1)+dydksi*f(i,j,2))*rho(i,j)/hksi(i,j)*a1

            dxdksi = 0.5D+00*(xc(i+2,j+1)+xc(i+2,j)-xc(i+1,j+1)-xc(i+1,j))
            dydksi = 0.5D+00*(yc(i+2,j+1)+yc(i+2,j)-yc(i+1,j+1)-yc(i+1,j))
            ruki1j = (dxdksi*f(i+1,j,1)+dydksi*f(i+1,j,2))*rho(i+1,j) &
              /hksi(i+1,j)*a2

            t1 = rukij*(apr*hksi(i+1,j)-ap(i,j)*hksi(i,j)/a1)
            t2 = ruki1j*(apr*hksi(i,j)-ap(i+1,j)*hksi(i+1,j)/a2)
            t3 = 0.5D+00*(t1+t2)

            if ( i == 2 ) then
              rukt = ruksi(2,j)
            else
              rukt = temp
            end if

            temp = ruksi(i+1,j)
            frc = hksi(i+1,j)/(hksi(i,j)+hksi(i+1,j))
            t4 = 0.5D+00*(frc*(ruksi(i+1,j)*ak1(i+1,j)-rukt*ak1(i,j)) &
              +(1.0D+00-frc)*(ruksi(i+1,j)*ak1(i+1,j) &
              -ruksi(i+2,j)*ak1(i+2,j)))

            ruksi(i+1,j) = ((hksip(i,j)+hksip(i+1,j)+t3)/denom+t4)/ak1(i+1,j)

            aip(i,j) = 1.0D+00/denom
            aim(i+1,j) = aip(i,j)

          end do

          aip(l1-1,j) = 0.0D+00
          aim(2,j) = 0.0D+00

        end do

        do i = 2,l1-1
          do j = 2,m1-2

            a1 = 0.5D+00*(ae1(i,j)+ae1(i,j+1))
            a2 = 0.5D+00*(ae1(i,j+1)+ae1(i,j+2))
            denom = 0.5D+00*(ap(i,j)*heta(i,j)/a1+ap(i,j+1)*heta(i,j+1)/a2)
            apr = 2.0D+00*denom/(heta(i,j)+heta(i,j+1))

            dxdeta = 0.5D+00*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5D+00*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            rueij = (dxdeta*f(i,j,1)+dydeta*f(i,j,2))*rho(i,j)/heta(i,j)*a1

            dxdeta = 0.5D+00*(xc(i+1,j+2)+xc(i,j+2)-xc(i+1,j+1)-xc(i,j+1))
            dydeta = 0.5D+00*(yc(i+1,j+2)+yc(i,j+2)-yc(i+1,j+1)-yc(i,j+1))
            rueij1 = (dxdeta*f(i,j+1,1)+dydeta*f(i,j+1,2))*rho(i,j+1)/heta(i,j+1)*a2

            t1 = rueij*(apr*heta(i,j+1)-ap(i,j)*heta(i,j)/a1)
            t2 = rueij1*(apr*heta(i,j)-ap(i,j+1)*heta(i,j+1)/a2)
            t3 = 0.5D+00*(t1+t2)
            ruet = rueta(i,2)
            if ( j /= 2) ruet = temp
            temp = rueta(i,j+1)
            frc = heta(i,j+1)/(heta(i,j)+heta(i,j+1))
            t4 = 0.5D+00*(frc*(rueta(i,j+1)*ae1(i,j+1)-ruet*ae1(i,j)) &
              +(1.0D+00-frc)*(rueta(i,j+1)*ae1(i,j+1) &
              -rueta(i,j+2)*ae1(i,j+2)))

            rueta(i,j+1) = ((hetap(i,j)+hetap(i,j+1)+t3)/denom+t4)/ae1(i,j+1)

            ajp(i,j) = 1.0D+00/denom
            ajm(i,j+1) = ajp(i,j)

          end do

          ajp(i,m1-1) = 0.0D+00
          ajm(i,2) = 0.0D+00

        end do

        do i = 2, l1-1
          do j = 2, m1-1
            con(i,j) = con(i,j)+ak1(i,j)*ruksi(i,j) &
              -ak1(i+1,j)*ruksi(i+1,j) &
              +ae1(i,j)*rueta(i,j)-ae1(i,j+1)*rueta(i,j+1)
            ap(i,j) = aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        smax = 0.0D+00
        ssum = 0.0D+00
        do i = 2,l1-1
          do j = 2,m1-1
            soor = ap(i,j)*f(i,j,3)-aip(i,j)*f(i+1,j,3) &
              -aim(i,j)*f(i-1,j,3)-ajp(i,j)*f(i,j+1,3) &
              -ajm(i,j)*f(i,j-1,3)-con(i,j)
            ssum = ssum+soor/rho(i,j)
            smax = max(abs(soor)/1.0D+04,smax)
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-1
            ap(i,j) = ap(i,j)/relax(np)
            con(i,j) = con(i,j)+(1.0D+00-relax(np))*ap(i,j)*f(i,j,3)
          end do
        end do

        call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
          heta,hksi,icrys,izone,jcrys,l1,m1,nf)

        call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)

        do j = 2,m1-1
          do i = 2,l1-2
            ruksi(i+1,j) = ruksi(i+1,j)+aip(i,j)/ak1(i+1,j) &
              *(f(i,j,nf)-f(i+1,j,nf))
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-2
            rueta(i,j+1) = rueta(i,j+1)+ajp(i,j)/ae1(i,j+1) &
              *(f(i,j,nf)-f(i,j+1,nf))
          end do
        end do

        do j = 2,m1-1

          do i = 2,l1-2

            bpi = f(i,j,3)
            bpi1 = f(i+1,j,3)-hksip(i,j)-hksip(i+1,j)
            em = hksi(i,j)*tmp(i,j)/(ak1(i,j)+ak1(i+1,j))
            ep = hksi(i+1,j)*tmp(i+1,j)/(ak1(i+1,j)+ak1(i+2,j))

            t1 = 0.5D+00*em*ep*(ruksi(i+2,j)*ak1(i+2,j) &
              -ruksi(i,j)*ak1(i,j))/(em+ep)
            bpm = (bpi*ep+bpi1*em)/(em+ep)+t1
            pksi(i+1,j) = bpm+hksip(i,j)

          end do

          pksi(l1,j) = 2.0D+00*f(l1-1,j,nf)-pksi(l1-1,j)
          pksi(2,j) = 2.0D+00*f(2,j,nf)-pksi(3,j)
          f(1,j,nf) = pksi(2,j)
          f(l1,j,nf) = pksi(l1,j)

        end do

        do i = 2,l1-1
          do j = 2,m1-2
            bpj = f(i,j,3)
            bpj1 = f(i,j+1,3)-hetap(i,j)-hetap(i,j+1)
            em = heta(i,j)*tmp(i,j)/(ae1(i,j)+ae1(i,j+1))
            ep = heta(i,j+1)*tmp(i,j+1)/(ae1(i,j+1)+ae1(i,j+2))

            t1 = 0.5D+00*em*ep*(rueta(i,j+2)*ae1(i,j+2) &
              -rueta(i,j)*ae1(i,j))/(em+ep)
            bpm = (bpj*ep+bpj1*em)/(em+ep)+t1
            peta(i,j+1) = bpm+hetap(i,j)

          end do

          peta(i,m1) = 2.0D+00 *f(i,m1-1,nf)-peta(i,m1-1)
          peta(i,2) = 2.0D+00 *f(i,2,nf)-peta(i,3)
          f(i,1,nf) = peta(i,2)
          f(i,m1,nf) = peta(i,m1)

        end do

        do j = 1,m1
          gam(1,j) = biu(j)
          gam(l1,j) = blu(j)
        end do

        do i = 1,l1
          gam(i,1) = bju(i)
          gam(i,m1) = bmu(i)
        end do

        do i = 1,l1
          do j = 1,m1
            con(i,j) = conu(i,j)
            ap(i,j) = cofu(i,j,1)
            aip(i,j) = cofu(i,j,2)
            aim(i,j) = cofu(i,j,3)
            ajp(i,j) = cofu(i,j,4)
            ajm(i,j) = cofu(i,j,5)
          end do
        end do

        do i = 2,l1-1

          if ( gam(i,1) == 0.0D+00 ) then
            ajm(i,2) = 0.0D+00
          end if

          if ( gam(i,m1) == 0.0D+00 ) then
            ajp(i,m1-1) = 0.0D+00
          end if

        end do

        do j = 2,m1-1

          if ( gam(1,j) == 0.0D+00 ) then
            aim(2,j) = 0.0D+00
          end if

          if ( gam(l1,j) == 0.0D+00 ) then
            aip(l1-1,j) = 0.0D+00
          end if

        end do

        do i = 2,l1-1
          do j = 2,m1-1

            xr = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j))
            yr = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j))
            xl = 0.5D+00*(xc(i,j+1)+xc(i,j))
            yl = 0.5D+00*(yc(i,j+1)+yc(i,j))
            xu = 0.5D+00*(xc(i+1,j+1)+xc(i,j+1))
            yu = 0.5D+00*(yc(i+1,j+1)+yc(i,j+1))
            xd = 0.5D+00*(xc(i+1,j)+xc(i,j))
            yd = 0.5D+00*(yc(i+1,j)+yc(i,j))
            dxdksi = xr-xl
            dydksi = yr-yl
            dxdeta = xu-xd
            dydeta = yu-yd
            dpksi = pksi(i+1,j)-pksi(i,j)
            dpeta = peta(i,j+1)-peta(i,j)
            t1 = -dydeta*dpksi+dydksi*dpeta
            t2 = dxdeta*dpksi-dxdksi*dpeta

            if ( mode == 1 ) then
              t1 = t1*r(i,j)
              t2 = t2*r(i,j)
            end if

            tmp(i,j) = t2
            con(i,j) = con(i,j)*vol(i,j)+t1
            ap(i,j) = -ap(i,j)*vol(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)

          end do
        end do

        nf = 1
        isol = 1

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          b3jbl,cappa,cd,cmu, &
          con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        do j = 1,m1
          gam(1,j) = biv(j)
          gam(l1,j) = blv(j)
        end do

        do i = 1,l1
          gam(i,1) = bjv(i)
          gam(i,m1) = bmv(i)
        end do

        do i = 1,l1
          do j = 1,m1
            con(i,j) = conv(i,j)
            ap(i,j) = cofv(i,j,1)
            aip(i,j) = cofv(i,j,2)
            aim(i,j) = cofv(i,j,3)
            ajp(i,j) = cofv(i,j,4)
            ajm(i,j) = cofv(i,j,5)
          end do
        end do

        nf = 2

        do i = 2,l1-1

          if ( gam(i,1) == 0.0D+00 ) then
            ajm(i,2) = 0.0D+00
          end if

          if ( gam(i,m1) == 0.0D+00 ) then
            ajp(i,m1-1) = 0.0D+00
          end if

        end do

        do j = 2,m1-1

          if ( gam(1,j) == 0.0D+00 ) then
            aim(2,j) = 0.0D+00
          end if

          if ( gam(l1,j) == 0.0D+00 ) then
            aip(l1-1,j) = 0.0D+00
          end if

        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = con(i,j)*vol(i,j)+tmp(i,j)
            ap(i,j) = -ap(i,j)*vol(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        isol = 1

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          b3jbl,cappa,cd,cmu, &
          con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)
!
!  N = 4, PC CORRECTED PRESSURE.
!
      else if ( n == 4 ) then

        do j = 1,m1
          gam(1,j) = biu(j)
          gam(l1,j) = blu(j)
        end do

        do i = 1,l1
          gam(i,1) = bju(i)
          gam(i,m1) = bmu(i)
        end do

        do i = 1,l1
          do j = 1,m1
            con(i,j) = conu(i,j)
            ap(i,j) = cofu(i,j,1)
            aip(i,j) = cofu(i,j,2)
            aim(i,j) = cofu(i,j,3)
            ajp(i,j) = cofu(i,j,4)
            ajm(i,j) = cofu(i,j,5)
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,1))*vol(i,j)
          end do
        end do

        nf = 1

        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          b3jbl,cappa,cd,cmu, &
          con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        nf = npc

        do i = 2,l1-1
          do j = 2,m1-1
            hutop(i,j) = (aip(i,j)*f(i+1,j,1)+aim(i,j)*f(i-1,j,1) &
              +ajp(i,j)*f(i,j+1,1)+ajm(i,j)*f(i,j-1,1)+con(i,j))/vol(i,j)
          end do
        end do

        do j = 2,m1-1
          do i = 2,l1-1
            con(i,j) = 0.0D+00
            ap(i,j) = 0.0D+00
          end do
        end do

        do j = 1,m1
          gam(1,j) = biv(j)
          gam(l1,j) = blv(j)
        end do

        do i = 1,l1
          gam(i,1) = bjv(i)
          gam(i,m1) = bmv(i)
        end do

        do i = 1,l1
          do j = 1,m1
            con(i,j) = conv(i,j)
            ap(i,j) = cofv(i,j,1)
            aip(i,j) = cofv(i,j,2)
            aim(i,j) = cofv(i,j,3)
            ajp(i,j) = cofv(i,j,4)
            ajm(i,j) = cofv(i,j,5)
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = (con(i,j)+ap(i,j)*f(i,j,2))*vol(i,j)
          end do
        end do

        nf = 2
        isol = 0

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          b3jbl,cappa,cd,cmu, &
          con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

        nf = npc

        do i = 2,l1-1
          do j = 2,m1-1
            hvtop(i,j) = (aip(i,j)*f(i+1,j,2)+aim(i,j)*f(i-1,j,2) &
              +ajp(i,j)*f(i,j+1,2)+ajm(i,j)*f(i,j-1,2) &
              +con(i,j))/vol(i,j)
            dxdksi = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            hksip(i,j) = 0.5D+00*(dxdksi*hutop(i,j)+dydksi*hvtop(i,j))
            dxdeta = 0.5D+00*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5D+00*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            hetap(i,j) = 0.5D+00*(dxdeta*hutop(i,j)+dydeta*hvtop(i,j))
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = 0.0D+00
            ap(i,j) = (aip(i,j)+aim(i,j)+ajm(i,j)+ajp(i,j))/vol(i,j)/rho(i,j)
            tmp(i,j) = ap(i,j)
          end do
        end do

        if ( .not.lortho ) then

          do j = 2,m1-1

            do i = 1,l1
              qt(i) = 0.5D+00*(rueta(i,j)+rueta(i,j+1))
            end do

            do i = 2,l1-1
              tem1 = hksi(i+1,j)+hksi(i,j)
              tem2 = hksi(i-1,j)+hksi(i,j)
              t1 = hksi(i,j)/tem1
              t2 = hksi(i+1,j)/tem1
              t3 = hksi(i,j)/tem2
              t4 = hksi(i-1,j)/tem2
              con(i,j) = con(i,j)+(t1*qt(i+1)+t2*qt(i))*ak2(i+1,j) &
                -(t4*qt(i)+t3*qt(i-1))*ak2(i,j)
            end do

          end do

          do i = 2,l1-1

            do j = 1,m1
              qt(j) = 0.5D+00*(ruksi(i,j)+ruksi(i+1,j))
            end do

            do j = 2,m1-1
              tem1 = heta(i,j+1)+heta(i,j)
              tem2 = heta(i,j-1)+heta(i,j)
              t1 = heta(i,j)/tem1
              t2 = heta(i,j+1)/tem1
              t3 = heta(i,j)/tem2
              t4 = heta(i,j-1)/tem2
              con(i,j) = con(i,j)+(t1*qt(j+1)+t2*qt(j))*ae2(i,j+1) &
                -(t4*qt(j)+t3*qt(j-1))*ae2(i,j)
            end do
          end do

        end if

        do j = 2,m1-1

          temp = 0.0D+00

          do i = 2,l1-2

            a1 = 0.5D+00*(ak1(i,j)+ak1(i+1,j))
            a2 = 0.5D+00*(ak1(i+1,j)+ak1(i+2,j))
            denom = 0.5D+00*(ap(i,j)*hksi(i,j)/a1+ap(i+1,j)*hksi(i+1,j)/a2)
            apr = 2.0D+00*denom/(hksi(i,j)+hksi(i+1,j))

            dxdksi = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j)-xc(i,j+1)-xc(i,j))
            dydksi = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j)-yc(i,j+1)-yc(i,j))
            rukij = (dxdksi*f(i,j,1)+dydksi*f(i,j,2))*rho(i,j)/hksi(i,j)*a1

            dxdksi = 0.5D+00*(xc(i+2,j+1)+xc(i+2,j)-xc(i+1,j+1)-xc(i+1,j))
            dydksi = 0.5D+00*(yc(i+2,j+1)+yc(i+2,j)-yc(i+1,j+1)-yc(i+1,j))
            ruki1j = (dxdksi*f(i+1,j,1)+dydksi*f(i+1,j,2))*rho(i+1,j) &
              /hksi(i+1,j)*a2

            t1 = rukij*(apr*hksi(i+1,j)-ap(i,j)*hksi(i,j)/a1)
            t2 = ruki1j*(apr*hksi(i,j)-ap(i+1,j)*hksi(i+1,j)/a2)
            t3 = 0.5D+00*(t1+t2)

            if ( i == 2 ) then
              rukt = ruksi(2,j)
            else
              rukt = temp
            end if

            temp = ruksi(i+1,j)
            frc = hksi(i+1,j)/(hksi(i,j)+hksi(i+1,j))
            t4 = 0.5D+00*(frc*(ruksi(i+1,j)*ak1(i+1,j)-rukt*ak1(i,j)) &
              +(1.0D+00-frc)*(ruksi(i+1,j)*ak1(i+1,j) &
              -ruksi(i+2,j)*ak1(i+2,j)))

            ruksi(i+1,j) = ((hksip(i,j)+hksip(i+1,j)+t3)/denom+t4)/ak1(i+1,j)

            ruksi(i+1,j) = ruksi(i+1,j)+(f(i,j,3)-f(i+1,j,3))/denom/ak1(i+1,j)

            aip(i,j) = 1.0D+00/denom
            aim(i+1,j) = aip(i,j)

          end do

          aip(l1-1,j) = 0.0D+00
          aim(2,j) = 0.0D+00

        end do

        do i = 2,l1-1
          do j = 2,m1-2

            a1 = 0.5D+00*(ae1(i,j)+ae1(i,j+1))
            a2 = 0.5D+00*(ae1(i,j+1)+ae1(i,j+2))
            denom = 0.5D+00*(ap(i,j)*heta(i,j)/a1+ap(i,j+1)*heta(i,j+1)/a2)
            apr = 2.0D+00*denom/(heta(i,j)+heta(i,j+1))

            dxdeta = 0.5D+00*(xc(i+1,j+1)+xc(i,j+1)-xc(i+1,j)-xc(i,j))
            dydeta = 0.5D+00*(yc(i+1,j+1)+yc(i,j+1)-yc(i+1,j)-yc(i,j))
            rueij = (dxdeta*f(i,j,1)+dydeta*f(i,j,2))*rho(i,j)/heta(i,j)*a1

            dxdeta = 0.5D+00*(xc(i+1,j+2)+xc(i,j+2)-xc(i+1,j+1)-xc(i,j+1))
            dydeta = 0.5D+00*(yc(i+1,j+2)+yc(i,j+2)-yc(i+1,j+1)-yc(i,j+1))
            rueij1 = (dxdeta*f(i,j+1,1)+dydeta*f(i,j+1,2))*rho(i,j+1)/heta(i,j+1)*a2

            t1 = rueij*(apr*heta(i,j+1)-ap(i,j)*heta(i,j)/a1)
            t2 = rueij1*(apr*heta(i,j)-ap(i,j+1)*heta(i,j+1)/a2)
            t3 = 0.5D+00*(t1+t2)
            ruet = rueta(i,2)
            if ( j /= 2) ruet = temp
            temp = rueta(i,j+1)
            frc = heta(i,j+1)/(heta(i,j)+heta(i,j+1))
            t4 = 0.5D+00*(frc*(rueta(i,j+1)*ae1(i,j+1)-ruet*ae1(i,j)) &
              +(1.0D+00-frc)*(rueta(i,j+1)*ae1(i,j+1) &
              -rueta(i,j+2)*ae1(i,j+2)))

            rueta(i,j+1) = ((hetap(i,j)+hetap(i,j+1)+t3)/denom+t4)/ae1(i,j+1)

            rueta(i,j+1) = rueta(i,j+1) &
              +(f(i,j,3)-f(i,j+1,3))/denom/ae1(i,j+1)

            ajp(i,j) = 1.0D+00/denom
            ajm(i,j+1) = ajp(i,j)

          end do

          ajp(i,m1-1) = 0.0D+00
          ajm(i,2) = 0.0D+00

        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = con(i,j)+ak1(i,j)*ruksi(i,j)-ak1(i+1,j)*ruksi(i+1,j) &
              +ae1(i,j)*rueta(i,j)-ae1(i,j+1)*rueta(i,j+1)
            ap(i,j) = aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        call solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
          heta,hksi,icrys,izone,jcrys,l1,m1,nf)

        call solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)

        do j = 2,m1-1
          do i = 2,l1-2
            ruksi(i+1,j) = ruksi(i+1,j)+aip(i,j)/ak1(i+1,j) &
              *(f(i,j,nf)-f(i+1,j,nf))
          end do
        end do

        do i = 2,l1-1
          do j = 2,m1-2
            rueta(i,j+1) = rueta(i,j+1)+ajp(i,j)/ae1(i,j+1) &
              *(f(i,j,nf)-f(i,j+1,nf))
          end do
        end do

        do j = 2,m1-1

          do i = 2,l1-2

            bpi = f(i,j,3)
            bpi1 = f(i+1,j,3)-hksip(i,j)-hksip(i+1,j)
            em = hksi(i,j)*tmp(i,j)/(ak1(i,j)+ak1(i+1,j))
            ep = hksi(i+1,j)*tmp(i+1,j)/(ak1(i+1,j)+ak1(i+2,j))
            pksi(i+1,j) = (f(i,j,nf)*ep+f(i+1,j,nf)*em)/(em+ep)

          end do

          pksi(l1,j) = 2.0D+00*f(l1-1,j,nf)-pksi(l1-1,j)
          pksi(2,j) = 2.0D+00*f(2,j,nf)-pksi(3,j)
          f(1,j,nf) = pksi(2,j)
          f(l1,j,nf) = pksi(l1,j)

        end do

        do i = 2,l1-1
          do j = 2,m1-2

            bpj = f(i,j,3)
            bpj1 = f(i,j+1,3)-hetap(i,j)-hetap(i,j+1)
            em = heta(i,j)*tmp(i,j)/(ae1(i,j)+ae1(i,j+1))
            ep = heta(i,j+1)*tmp(i,j+1)/(ae1(i,j+1)+ae1(i,j+2))
            peta(i,j+1) = (f(i,j,nf)*ep+f(i,j+1,nf)*em)/(em+ep)

          end do

          peta(i,m1) = 2.0D+00*f(i,m1-1,nf)-peta(i,m1-1)
          peta(i,2) = 2.0D+00*f(i,2,nf)-peta(i,3)
          f(i,1,nf) = peta(i,2)
          f(i,m1,nf) = peta(i,m1)

        end do

        do i = 2,l1-1
          do j = 2,m1-1

            xr = 0.5D+00*(xc(i+1,j+1)+xc(i+1,j))
            yr = 0.5D+00*(yc(i+1,j+1)+yc(i+1,j))
            xl = 0.5D+00*(xc(i,j+1)+xc(i,j))
            yl = 0.5D+00*(yc(i,j+1)+yc(i,j))
            xu = 0.5D+00*(xc(i+1,j+1)+xc(i,j+1))
            yu = 0.5D+00*(yc(i+1,j+1)+yc(i,j+1))
            xd = 0.5D+00*(xc(i+1,j)+xc(i,j))
            yd = 0.5D+00*(yc(i+1,j)+yc(i,j))
            dxdksi = xr-xl
            dydksi = yr-yl
            dxdeta = xu-xd
            dydeta = yu-yd
            dpksi = pksi(i+1,j)-pksi(i,j)
            dpeta = peta(i,j+1)-peta(i,j)
            t1 = -dydeta*dpksi+dydksi*dpeta
            t2 = dxdeta*dpksi-dxdksi*dpeta

            if ( mode == 1 ) then
              t1 = t1*r(i,j)
              t2 = t2*r(i,j)
            end if

            f(i,j,1) = f(i,j,1)+t1/tmp(i,j)/rho(i,j)/vol(i,j)
            f(i,j,2) = f(i,j,2)+t2/tmp(i,j)/rho(i,j)/vol(i,j)

          end do
        end do
!
!  N > 4, T, W, TK, TE, E, PSI
!
      else
!
!  Get the transport coefficients, GAM.
!
        call gamsor(ap,birad,ce1,ce2,cmu,con,ewall,f,fcsl,fksl, &
          fma,fo,frsl,gam,gamt,grash,hamag,heta,hksi,icrys,inturb, &
          izone,jcrys,l1,lsolve,m1,mode,nf,pr,r,rdtm,re,recb,rect, &
          rho,rhocon,rpr,sige,sigk,sigt,stel,stes,tal,tas,tf,tw,x, &
          xc,y,yc)

        do j = 2,m1-1

          diff = 2.0D+00*gam(1,j)/hksi(2,j)
          flow = -ruksi(2,j)
          call diflow(acof,diff,flow)
          aim(2,j) = acof*ak1(2,j)

          diff = 2.0D+00*gam(l1,j)/hksi(l1-1,j)
          flow = ruksi(l1,j)
          call diflow(acof,diff,flow)
          aip(l1-1,j) = acof*ak1(l1,j)

          do i = 2,l1-2
            diff = gam(i,j)*gam(i+1,j) &
              /(0.5D+00*hksi(i,j)*gam(i+1,j)+0.5D+00*hksi(i+1,j)*gam(i,j))
            flow = ruksi(i+1,j)
            call diflow(acof,diff,flow)
            aip(i,j) = acof*ak1(i+1,j)
            aim(i+1,j) = (acof+ruksi(i+1,j))*ak1(i+1,j)
          end do

        end do

        do i = 2,l1-1

          diff = 2.0D+00*gam(i,1)/heta(i,2)
          flow = -rueta(i,2)
          call diflow(acof,diff,flow)
          ajm(i,2) = acof*ae1(i,2)

          diff = 2.0D+00*gam(i,m1)/heta(i,m1-1)
          flow = rueta(i,m1)
          call diflow(acof,diff,flow)
          ajp(i,m1-1) = acof*ae1(i,m1)

          do j = 2,m1-2
            diff = gam(i,j)*gam(i,j+1) &
              /(0.5D+00*heta(i,j)*gam(i,j+1)+0.5D+00*heta(i,j+1)*gam(i,j))
            flow = rueta(i,j+1)
            call diflow(acof,diff,flow)
            ajp(i,j) = acof*ae1(i,j+1)
            ajm(i,j+1) = (acof+rueta(i,j+1))*ae1(i,j+1)
          end do

        end do

        do i = 2,l1-1
          do j = 2,m1-1
            con(i,j) = con(i,j)*vol(i,j)
            ap(i,j) = -ap(i,j)*vol(i,j)+aip(i,j)+aim(i,j)+ajp(i,j)+ajm(i,j)
          end do
        end do

        isol = 1

        call flux(ae1,ae2,aim,aip,ajm,ajp,ak1,ak2,ap,b1jbl,b2jbl, &
          b3jbl,cappa,cd,cmu, &
          con,epsil,f,fjeta,fjksi,fmax,fn,fo,gam,heta,hksi, &
          icrys,isol,iter,izone,jcrys,l1,lblk,lconv,lortho, &
          m1,mode,nf,nsolve,ntimes,r,relax,res,rueta,ruksi)

      end if

    end if

  end do

  return
end
subroutine setx ( l, m, ni, nj, x, xc, y, yc )

!*****************************************************************************80
!
!! SETX calculates the locations of the primary nodes.
!
!  Discussion:
!
!    SETX is given (XC,YC), the locations of the "corners" of the
!    control volumes, and calculates the locations of the primary nodes
!    (X,Y).
!
!    (XC,YC) sits in the center of the control volume.  On the other
!    hand, the point (X,Y) is NOT necessarily the center of the four
!    nearby values of XC and YC.
!
!    SETX may be called to do a portion of the region, or the
!    full region (in which case L = L0, M=M0).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, integer ( kind = 4 ) L.
!    L specifies the portion of the vertical grid size to use.
!
!    Input, integer ( kind = 4 ) M.
!    M specifies the portion of the horizontal grid size to use.
!
!    Input, integer ( kind = 4 ) NI.
!    NI is the maximum number of grid points in the "I" or vertical
!    direction.
!
!    Input, integer ( kind = 4 ) NJ.
!    NJ is the maximum number of grid points in the "J" or horizontal
!    direction.
!
!    Output, real ( kind = 8 ) X(NI,NJ).
!    X coordinate of primary nodes, which are the centers of
!    control volumes.
!    In most cases, X(I,J) is the average of the four enclosing
!    corner values, XC(I,J), XC(I+1,J), XC(I,J+1), and XC(I+1,J+1).
!
!    Input, real ( kind = 8 ) XC(NI,NJ).
!    XC contains the X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Output, real ( kind = 8 ) Y(NI,NJ).
!    Y coordinate of primary nodes, which are the centers of
!    control volumes.
!    In most cases, Y(I,J) is the average of the four enclosing
!    corner values, YC(I,J), YC(I+1,J), YC(I,J+1), and YC(I+1,J+1).
!
!    Input, real ( kind = 8 ) YC(NI,NJ).
!    YC contains the Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ) ni
  integer ( kind = 4 ) nj

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) y(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
!
!  The central set of values.
!
  do i = 2,l-1
    do j = 2,m-1
      x(i,j) = 0.25D+00*(xc(i,j)+xc(i+1,j)+xc(i,j+1)+xc(i+1,j+1))
      y(i,j) = 0.25D+00*(yc(i,j)+yc(i+1,j)+yc(i,j+1)+yc(i+1,j+1))
    end do
  end do
!
!  The first and last columns, I = 2 to L-1, J=1 and J=M.
!
  do i = 2,l-1
    x(i,1) = 0.5D+00*(xc(i,2)+xc(i+1,2))
    y(i,1) = 0.5D+00*(yc(i,2)+yc(i+1,2))
    x(i,m) = 0.5D+00*(xc(i,m)+xc(i+1,m))
    y(i,m) = 0.5D+00*(yc(i,m)+yc(i+1,m))
  end do
!
!  The first and last rows, I = 1 and I=L, J=2 to M-1.
!
  do j = 2,m-1
    x(1,j) = 0.5D+00*(xc(2,j)+xc(2,j+1))
    y(1,j) = 0.5D+00*(yc(2,j)+yc(2,j+1))
    x(l,j) = 0.5D+00*(xc(l,j)+xc(l,j+1))
    y(l,j) = 0.5D+00*(yc(l,j)+yc(l,j+1))
  end do
!
!  The corners.
!
  x(1,1) = xc(2,2)
  y(1,1) = yc(2,2)

  x(1,m) = xc(2,m)
  y(1,m) = yc(2,m)

  x(l,1) = xc(l,2)
  y(l,1) = yc(l,2)

  x(l,m) = xc(l,m)
  y(l,m) = yc(l,m)

  return
end
subroutine solve1(aim,aip,ajm,ajp,ap,cappa,cd,cmu,con,f,fo, &
  heta,hksi,icrys,izone,jcrys,l1,m1,nf)

!*****************************************************************************80
!
!! SOLVE1 sets or modifies certain entries of the linear system.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: ns = 10

  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) cappa
  real ( kind = 8 ) cd
  real ( kind = 8 ) cm4
  real ( kind = 8 ) cmu
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) f(ni,nj,ns)
  real ( kind = 8 ) fo(ni,nj,ns)
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) izone
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcrys
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) nf
!
!  Temperature calculation in solid zone.
!
  if ( izone == 1 .and. nf == 5 ) then
    do i = 1,icrys
      do j = 1,m1
        ap(i,j) = 1.0D+00
        aip(i,j) = 0.0D+00
        aim(i,j) = 0.0D+00
        ajp(i,j) = 0.0D+00
        ajm(i,j) = 0.0D+00
        con(i,j) = fo(i,j,5)
      end do
    end do
  end if
!
!  Magnetic stream function calculation in solid zone.
!
  if ( izone == 1 .and. nf == 9 ) then
    do i = 1,icrys-1
      do j = 1,m1
        ap(i,j) = 1.0D+00
        aip(i,j) = 0.0D+00
        aim(i,j) = 0.0D+00
        ajp(i,j) = 0.0D+00
        ajm(i,j) = 0.0D+00
        con(i,j) = fo(i,j,9)
      end do
    end do
  end if
!
!  U, V, P, PC, or T calculation in liquid zone.
!
  if ( izone == 2 .and. nf <= 5 ) then
    do j = 2,jcrys
      ap(icrys,j) = 1.0D+00
      aip(icrys,j) = 0.0D+00
      aim(icrys,j) = 0.0D+00
      ajp(icrys,j) = 0.0D+00
      ajm(icrys,j) = 0.0D+00
      con(icrys,j) = 0.0D+00
    end do
  end if
!
!  Turbulent dissipation (TE) calculation.
!
  if ( nf == 8 ) then

    cm4 = cmu**0.25

    do i = 2,l1-1
      con(i,2) = cd*f(i,2,7)**1.5/(cappa*cm4*0.5D+00*heta(i,2))
      ap(i,2) = 1.0D+00
      con(i,m1-1) = cd*f(i,m1-1,7)**1.5/(cappa*cm4*0.5D+00*heta(i,m1-1))
      ap(i,m1-1) = 1.0D+00
    end do

    do j = 2,m1-1
      con(2,j) = cd*f(2,j,7)**1.5/(cappa*cm4*0.5D+00*hksi(2,j))
      ap(2,j) = 1.0D+00
      con(l1-1,j) = cd*f(l1-1,j,7)**1.5/(cappa*cm4*0.5D+00*hksi(l1-1,j))
      ap(l1-1,j) = 1.0D+00
    end do

  end if

  return
end
subroutine solve2(aim,aip,ajm,ajp,ap,con,f,l1,lblk,m1,nf,nsolve)

!*****************************************************************************80
!
!! SOLVE2 is the tridiagonal matrix solver.
!
!  Discussion:
!
!    SOLVE2 must solve a set of linear equations defined on a five point
!    stencil, of the form:
!
!      A(I,J)*  U(I,J)
!    - A(I-1,J)*U(I-1,J) - A(I+1,J)*U(I+1,J)
!    - A(I,J-1)*U(I,J-1) - A(I,J+1)*U(I,J+1)  =  CON.
!
!    The equations must be solved for indices 2 < =  I < L1-1
!    and 2 <= J <= M1-1.
!
!    If I = 1 or I=L1 or J=1 or J=M1, then U(I,J) already contains its
!    correct value, and these correct values must be included in the
!    above equations.
!
!    The solution of the equations is stored in F(*,*,NF).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64
  integer ( kind = 4 ), parameter :: nmaxij = 64
  integer ( kind = 4 ), parameter :: ns = 10

  real ( kind = 8 ) aim(ni,nj)
  real ( kind = 8 ) aip(ni,nj)
  real ( kind = 8 ) ajm(ni,nj)
  real ( kind = 8 ) ajp(ni,nj)
  real ( kind = 8 ) ap(ni,nj)
  real ( kind = 8 ) bl
  real ( kind = 8 ) blc
  real ( kind = 8 ) blm
  real ( kind = 8 ) blp
  real ( kind = 8 ) con(ni,nj)
  real ( kind = 8 ) denom
  real ( kind = 8 ) f(ni,nj,ns)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  logical lblk(ns)
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nsolve(ns)
  integer ( kind = 4 ) nt
  real ( kind = 8 ) pt(nmaxij)
  real ( kind = 8 ) qt(nmaxij)

  do nt = 1,nsolve(nf)

    if ( lblk(nf) ) then

      pt(1) = 0.0D+00
      qt(1) = 0.0D+00

      do i = 2,l1-1

        bl = 0.0D+00
        blp = 0.0D+00
        blm = 0.0D+00
        blc = 0.0D+00

        do j = 2,m1-1
          bl = bl+ap(i,j)
          if ( j /= m1-1) bl = bl-ajp(i,j)
          if ( j /= 2) bl = bl-ajm(i,j)
          blp = blp+aip(i,j)
          blm = blm+aim(i,j)
          blc = blc+con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
            +ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf)-ap(i,j)*f(i,j,nf)
        end do

        denom = bl-pt(i-1)*blm

        if ( abs(denom/bl) < 1.0D-10 ) then
          pt(i) = 0.0D+00
          qt(i) = 0.0D+00
        else
          pt(i) = blp/denom
          qt(i) = (blc+blm*qt(i-1))/denom
        end if

      end do

      bl = 0.0D+00
      do i = l1-1,2,-1

        bl = bl*pt(i)+qt(i)

        do j = 2,m1-1
          f(i,j,nf) = f(i,j,nf)+bl
        end do

      end do

      pt(1) = 0.0D+00
      qt(1) = 0.0D+00

      do j = 2,m1-1

        bl = 0.0D+00
        blp = 0.0D+00
        blm = 0.0D+00
        blc = 0.0D+00

        do i = 2,l1-1
          bl = bl+ap(i,j)
          if ( i /= l1-1) bl = bl-aip(i,j)
          if ( i /= 2) bl = bl-aim(i,j)
          blp = blp+ajp(i,j)
          blm = blm+ajm(i,j)
          blc = blc+con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
            +ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf)-ap(i,j)*f(i,j,nf)
        end do

        denom = bl-pt(j-1)*blm

        if ( abs(denom/bl) < 1.0D-10 ) then
          pt(j) = 0.0D+00
          qt(j) = 0.0D+00
        else
          pt(j) = blp/denom
          qt(j) = (blc+blm*qt(j-1))/denom
        end if

      end do

      bl = 0.0D+00
      do j = m1-1,2,-1
        bl = bl*pt(j)+qt(j)
        do i = 2,l1-1
          f(i,j,nf) = f(i,j,nf)+bl
        end do
      end do

    end if
!
!  Sweep 1: For each column J increasing.
!
    do j = 2,m1-1

      pt(1) = 0.0D+00
      qt(1) = f(1,j,nf)

      do i = 2,l1-1
        pt(i) = aip(i,j)/(ap(i,j)-pt(i-1)*aim(i,j))
        qt(i) = (con(i,j)+ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf) &
          +aim(i,j)*qt(i-1))/(ap(i,j)-pt(i-1)*aim(i,j))
      end do

      do i = l1-1,2,-1
        f(i,j,nf) = f(i+1,j,nf)*pt(i)+qt(i)
      end do

    end do
!
!  Sweep 2: For each column J decreasing.
!
    do j = m1-2,2,-1

      pt(1) = 0.0D+00
      qt(1) = f(1,j,nf)

      do i = 2,l1-1
        pt(i) = aip(i,j)/(ap(i,j)-pt(i-1)*aim(i,j))
        qt(i) = (con(i,j)+ajp(i,j)*f(i,j+1,nf)+ajm(i,j)*f(i,j-1,nf) &
          +aim(i,j)*qt(i-1))/(ap(i,j)-pt(i-1)*aim(i,j))
      end do

      do i = l1-1,2,-1
        f(i,j,nf) = f(i+1,j,nf)*pt(i)+qt(i)
      end do

    end do
!
!  Sweep 3: For each row I increasing.
!
    do i = 2,l1-1

      pt(1) = 0.0D+00
      qt(1) = f(i,1,nf)

      do j = 2,m1-1
        pt(j) = ajp(i,j)/(ap(i,j)-pt(j-1)*ajm(i,j))
        qt(j) = (con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
          +ajm(i,j)*qt(j-1))/(ap(i,j)-pt(j-1)*ajm(i,j))
      end do

      do j = m1-1,2,-1
        f(i,j,nf) = f(i,j+1,nf)*pt(j)+qt(j)
     end do

    end do
!
!  Sweep 4: For each row I decreasing.
!
    do i = l1-1,2,-1

      pt(1) = 0.0D+00
      qt(1) = f(i,1,nf)

      do j = 2,m1-1
        pt(j) = ajp(i,j)/(ap(i,j)-pt(j-1)*ajm(i,j))
        qt(j) = (con(i,j)+aip(i,j)*f(i+1,j,nf)+aim(i,j)*f(i-1,j,nf) &
          +ajm(i,j)*qt(j-1))/(ap(i,j)-pt(j-1)*ajm(i,j))
      end do

      do j = m1-1,2,-1
        f(i,j,nf) = f(i,j+1,nf)*pt(j)+qt(j)
      end do

    end do

  end do

  return
end
subroutine vortic ( heta, hksi, icrys, l0, m0, u, v, vort, xc, yc )

!*****************************************************************************80
!
!! VORTIC computes the value of the vorticity.
!
!  Discussion:
!
!    VORTIC estimates the value of the vorticity function VORT(I,J) at
!    each of the primary nodes (I,J), based on the values of the
!    horizontal and vertical velocities U and V, the XSI and ETA spacings
!    of the control volumes, and the X and Y coordinates of the
!    corners of the control volumes.
!
!    The fluid flow vorticity is defined as
!
!      VORT(X,Y)  =  dV(X,Y)/dX - dU(X,Y)/dY
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real ( kind = 8 ) HETA(NI,NJ).
!    HETA is the physical length of the control volume in the ETA direction.
!    This is the actual length of the line connecting the mid-side
!    nodes that separate primary node (I,J) from primary nodes
!    (I,J-1) and (I,J+1).
!
!    Input, real ( kind = 8 ) HKSI(NI,NJ).
!    HKSI is the physical length of the control volume in the KSI direction.
!    This is the actual length of the line connecting the mid-side
!    nodes that separate primary node (I,J) from primary nodes
!    (I-1,J) and (I+1,J).
!
!    Input, integer ( kind = 4 ) ICRYS.
!    ICRYS specifies the end of the crystal in the I array direction,
!    and in the vertical coordinate direction.
!
!    Input, integer ( kind = 4 ) L0.
!    L0 is the extent of the grid in the vertical or "I" coordinate.
!
!    Input, integer ( kind = 4 ) M0.
!    M0 is the extent of the grid in the horizontal or "J" coordinate.
!
!    Input, real ( kind = 8 ) U(NI,NJ).
!    U is the horizontal component of velocity at each primary node (I,J).
!
!    Input, real ( kind = 8 ) V(NI,NJ).
!    V is the vertical component of velocity at each primary node (I,J).
!
!    Output, real ( kind = 8 ) VORT(NI,NJ).
!    VORT is the fluid flow vorticity at each primary node (I,J).
!
!    Input, real ( kind = 8 ) XC(NI,NJ).
!    XC contains the X coordinate of nodes which are the "corners" of
!    control volumes.
!
!    Input, real ( kind = 8 ) YC(NI,NJ).
!    YC is the Y coordinate of nodes which are the "corners" of
!    control volumes.
!
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) heta(ni,nj)
  real ( kind = 8 ) hksi(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icrys
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) m0
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) vort(ni,nj)
  real ( kind = 8 ) xc(ni,nj)
  real ( kind = 8 ) yc(ni,nj)
!
!  At all the interior nodes, get the vorticitiy.
!  In the liquid region, get the spatial derivatives of the flow,
!  but in the void and crystal regions, just set it to zero.
!
  do i = 2, l0-1
    do j = 2, m0-1

      if ( i < icrys ) then
        call gradnt(heta,hksi,i,j,u,dudx,dudy,xc,yc)
        call gradnt(heta,hksi,i,j,v,dvdx,dvdy,xc,yc)
        vort(i,j) = dvdx-dudy
      else
        vort(i,j) = 0.0D+00
      end if

    end do
  end do
!
!  Copy the vorticity data to the sides and corners.
!
  do i = 2,l0-1
    vort(i,1) = vort(i,2)
    vort(i,m0) = vort(i,m0-1)
  end do

  do j = 2,m0
    vort(1,j) = vort(2,j)
    vort(l0,j) = vort(l0-1,j)
  end do

  vort(1,1) = vort(2,2)
  vort(1,m0) = vort(2,m0-1)
  vort(l0,1) = vort(l0-1,2)
  vort(l0,m0) = vort(l0-1,m0-1)

  return
end
subroutine wrtec ( gamt, l1, m1, p, psi, t, te, tk, u, v, w, x, y )

!*****************************************************************************80
!
!! WRTEC writes solution data to a TECPLOT input file.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
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
!    Input, real GAMT(NI,NJ).
!    GAMT(I,J) is the diffusion coefficient at primary node (I,J).
!
!    Input L1.
!    L1 is the maximum value for I when indexing the data arrays.
!
!    Input M1.
!    M1 is the maximum value for J when indexing the data arrays.
!
!    Input, real P(NI,NJ).
!    P(I,J) is the pressure at primary node (I,J).
!
!    Input, real PSI(NI,NJ).
!    PSI(I,J) is the stream function at primary node (I,J).
!
!    Input, real T(NI,NJ).
!    T(I,J) is the temperature at primary node (I,J).
!
!    Input, real TE(NI,NJ).
!    TE(I,J) is the turbulent epsilon coefficient at
!    primary node (I,J).
!
!    Input, real TK(NI,NJ).
!    TK(I,J) is the turbulent K coefficient at primary node (I,J).
!
!    Input, real U(NI,NJ).
!    U(I,J) is the horizontal velocity at primary node (I,J).
!
!    Input, real V(NI,NJ).
!    V(I,J) is the vertical velocity at primary node (I,J).
!
!    Input, real W(NI,NJ).
!    W(I,J) is the axial velocity at primary node (I,J).
!
!    Input, real X(NI,NJ).
!    X(I,J) is the X coordinate of primary node (I,J).
!
!    Input, real Y(NI,NJ).
!    Y(I,J) is the Y coordinate of primary node (I,J).
!
  implicit none

  integer ( kind = 4 ), parameter :: ni = 64
  integer ( kind = 4 ), parameter :: nj = 64

  character ( len = 11 ) filtec
  real ( kind = 8 ) gamt(ni,nj)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iunit
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) m1
  real ( kind = 8 ) p(ni,nj)
  real ( kind = 8 ) psi(ni,nj)
  real ( kind = 8 ) t(ni,nj)
  real ( kind = 8 ) te(ni,nj)
  real ( kind = 8 ) tk(ni,nj)
  real ( kind = 8 ) u(ni,nj)
  real ( kind = 8 ) v(ni,nj)
  real ( kind = 8 ) w(ni,nj)
  real ( kind = 8 ) x(ni,nj)
  real ( kind = 8 ) y(ni,nj)
!
!  Some local variables.
!
  filtec = 'tecplot.txt'
  iunit = 10
!
!  Delete any old copy of the TECPLOT data file.
!
  open(unit = iunit,file=filtec,status='old',err=10)
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRTEC - Note:'
  write ( *, '(a)' ) '  Deleting a previous copy of '//filtec
  close(unit = iunit,status='delete')

10    continue
!
!  Open a new copy of the TECPLOT data file.
!
  open ( unit = iunit, file = filtec, status = 'new' )
!
!  Write the data.
!  Note that the crystal program has X and Y reversed.
!  To keep TECPLOT from going mad, we reverse X and Y, and U and V.
!
  write ( iunit, '(a)' )' title = ','"zhanghui"'
  write ( iunit, '(a)' )' variables  = "x","y","u","v","w",'// &
    '"Press", "Temp","Psi","TK","TE","Gamt"'
  write ( iunit, '(a)' )' zone t =  "zone 1" , i=',l1,' , j=',m1,' , f= block'
  write ( iunit, '(5g15.6)' ) ((y(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((x(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((v(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((u(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((w(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((p(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((t(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((psi(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((tk(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((te(i,j),i = 1,l1),j=1,m1)
  write ( iunit, '(5g15.6)' ) ((gamt(i,j),i = 1,l1),j=1,m1)
!
!  Close the file.
!
  close ( unit = iunit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WRTEC - Note:'
  write ( *, '(a)' ) '  TECPLOT data written to "' // trim ( filtec ) // '".'

  return
end
