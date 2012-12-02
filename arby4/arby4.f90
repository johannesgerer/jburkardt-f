program main

!*****************************************************************************80
!
!! MAIN is the main program for ARBY4.
!
!  Discussion:
!
!    ARBY4 solves a fluid flow problem which has several parameters.
!
!    ARBY4 can solve the problem using the finite element method.  The
!    resulting solution will be called a "full" solution.
!
!    Once a full solution has been calculated, ARBY4 can compute the
!    first several derivatives of the solution with respect to the
!    Reynolds number, derive a reduced basis, and find reduced
!    solutions to the problem at a variety of Reynolds numbers.
!
!  Diary:
!
!    12 December 2000
!
!    Double precision constants REQUIRE D+00 or something similar, otherwise
!    they are computed as though they were REAL.  So I tried to fix that,
!    quickly.
!
!    14 September 1996
! 
!    Mr Lee wants a code that will solve the driven cavity problem.
!    I thought I would get FLOW4, but somehow ARBY4 may be the
!    quicker solution.
! 
!    12 September 1996
! 
!    Worked on GFL2RB and GRB2FL.
!  
!    11 September 1996
! 
!    Max urged me to write a reduced basis paper, as Kazi Ito has already
!    done so.
! 
!    I updated GFL2RB.
!  
!    17 August 1996
! 
!    OK, making transition BACK to representation of full solution as
!    GFL(GRB) = GFLRB + Q * GRB.
!  
!    16 August 1996
! 
!    GFLRB is already stored.  So now I have to reformulate the representation.
!  
!    15 August 1996
! 
!    I decided that I had to treat boundary conditions as follows:
! 
!      The reduced solution is represented as Y = Y0 + YBC*cbc + YFE*cfe
! 
!      Here, Y0 is the point at which the basis was generated.
!      Each column of YBC corresponds to the flow solution of the problem
!      at Y0 with boundary conditions differentiated with respect to
!      parameter I.  The YFE stuff is as usual.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none
!
!  Set parameters that are independent.
!
  integer ( kind = 4 ), parameter :: liv = 60
  integer ( kind = 4 ), parameter :: maxbcrb = 3
  integer ( kind = 4 ), parameter :: maxferb = 6
  integer ( kind = 4 ), parameter :: maxnx = 21
  integer ( kind = 4 ), parameter :: maxny = 21
  integer ( kind = 4 ), parameter :: maxparb = 5
  integer ( kind = 4 ), parameter :: maxparf = 5
!
!  Set parameters that are dependent on parameters.
!
!  The assignment of LDAFL should really read (ldafl = 29*min(nx,ny)).
!
  integer ( kind = 4 ), parameter :: ldafl = 29 * maxny
  integer ( kind = 4 ), parameter :: maxcofrb = maxbcrb + maxferb
  integer ( kind = 4 ), parameter :: maxelm = 2 * ( maxnx - 1 ) * ( maxny - 1 )
  integer ( kind = 4 ), parameter :: maxnfl = &
    2 * ( 2 * maxnx - 1 ) * ( 2 * maxny - 1 ) + maxnx * maxny
  integer ( kind = 4 ), parameter :: maxnp = ( 2 * maxnx - 1 ) * ( 2 * maxny - 1 )
  integer ( kind = 4 ), parameter :: maxpar = maxparb + maxparf + 1
!
!  Set parameters that are dependent on parameters that are dependent
!  on parameters.
!
  integer ( kind = 4 ), parameter :: lv = 78+maxpar*(maxpar+21)/2

  real ( kind = 8 ) afl(ldafl,maxnfl)
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,maxelm)
  character ( len = 9 ) chtime
  character ( len = 80 ) command
  real ( kind = 8 ) cost
  real ( kind = 8 ) cost0
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  character ( len = 8 ) date
  real ( kind = 8 ) detlog
  real ( kind = 8 ) detman
  real ( kind = 8 ) difcof(maxcofrb)
  character ( len = 30 ) disfil
  real ( kind = 8 ) dopt(maxpar)
  real ( kind = 8 ) dpar
  real ( kind = 8 ) drey
  logical dvneq
  logical echo
  real ( kind = 8 ) epsdif
  character ( len = 2 ) eqn(maxnfl)
  real estart
  real estop
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) factj
  real ( kind = 8 ) gfl(maxnfl)
  real ( kind = 8 ) gflafl(maxnfl)
  real ( kind = 8 ) gflnrm
  real ( kind = 8 ) gflopt(maxnfl)
  real ( kind = 8 ) gflrb(maxnfl)
  real ( kind = 8 ) gflsav(maxnfl)
  real ( kind = 8 ) gflsen(maxnfl)
  real ( kind = 8 ) gfltar(maxnfl)
  real ( kind = 8 ) gfltay(maxnfl)
  real ( kind = 8 ) gfltmp(maxnfl)
  real ( kind = 8 ) grb(maxcofrb)
  real ( kind = 8 ) grbarb(maxcofrb)
  real ( kind = 8 ) grbopt(maxcofrb)
  real ( kind = 8 ) grbsav(maxcofrb)
  real ( kind = 8 ) grbsen(maxcofrb)
  real ( kind = 8 ) grbtay(maxcofrb)
  character ( len = 20 ) gridx
  character ( len = 20 ) gridy
  real ( kind = 8 ) gsen(maxcofrb)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) icolrb(maxcofrb)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx(3,maxnp)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivfl(maxnfl)
  integer ( kind = 4 ) ipivrb(maxcofrb)
  integer ( kind = 4 ) iprint
  integer ( kind = 4 ) isotri(maxelm)
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ival1
  integer ( kind = 4 ) ival2
  integer ( kind = 4 ) ivopt(liv)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jtay
  integer ( kind = 4 ) klo
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) lenc
  logical s_eqi
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) maxsim
  integer ( kind = 4 ) mhi
  integer ( kind = 4 ) mlo
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) nlo
  integer ( kind = 4 ) node(6,maxelm)
  integer ( kind = 4 ) nodelm(maxnp)
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nprof(2*maxny-1)
  integer ( kind = 4 ) nsenfl
  integer ( kind = 4 ) ntay
  integer ( kind = 4 ) numdif
  integer ( kind = 4 ) numnew
  integer ( kind = 4 ) numopt
  integer ( kind = 4 ) numsim
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) p(maxnp)
  real ( kind = 8 ) par(maxpar)
  real ( kind = 8 ) parafl(maxpar)
  real ( kind = 8 ) pararb(maxpar)
  real ( kind = 8 ) pardif ( maxpar)
  real ( kind = 8 ) paropt(maxpar)
  real ( kind = 8 ) parrb(maxpar)
  real ( kind = 8 ) parsav(maxpar)
  real ( kind = 8 ) parsen(maxpar)
  real ( kind = 8 ) partar(maxpar)
  real ( kind = 8 ) phifl(3,6,10,maxelm)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(maxnfl)
  real ( kind = 8 ) resflsav(maxnfl)
  real ( kind = 8 ) resfltmp(maxnfl)
  real ( kind = 8 ) resrb(maxcofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) reytay
  real ( kind = 8 ) rbase(maxcofrb,maxcofrb)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) senrb(maxcofrb,maxcofrb)
  real ( kind = 8 ) splbmp(maxparb+2)
  real ( kind = 8 ) splflo(maxparf)
  real tarray(2)
  real ( kind = 8 ) taubmp(maxparb+2)
  real ( kind = 8 ) tauflo(maxparf)
  character ( len = 30 ) tecfil
  real ( kind = 8 ) temp
  character ( len = 10 ) time
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  character ( len = 10 ) tstart
  character ( len = 10 ) tstop
  real ( kind = 8 ) u(maxnp)
  real ( kind = 8 ) v(maxnp)
  real ( kind = 8 ) value
  real ( kind = 8 ) vopt(lv)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xopt(maxpar)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xquad(3,maxelm)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yquad(3,maxelm)
  real ( kind = 8 ) yrange

  call timestamp ( )

  iprint = 1
!
!  Get initial CPU clock reading.
!
  call etime ( tarray, estart )

  echo = .false.

  call hello ( maxnx, maxny )
!
!  Get the starting time.
!
  call date_and_time ( date, time )
  tstart = time
!
!  Open the file in which we record the user input.
!
  open ( unit = 17, file = 'arby.in', status = 'replace' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARBY4 - Init: Initialize all data.'

  call init(afl,arb,area,command,cost,costb,costp,costu, &
    costv,difcof,disfil,drey,epsdif,eqn,etaq,gfl,gflafl, &
    gflrb,gflsav,gflsen,gfltar,gfltay,grb,grbarb,grbsav, &
    grbsen,grbtay,gridx,gridy,hx,hy,ibs,ibump,icolrb,ierror, &
    ifs,ihi,ijac,ilo,indx,iopt,ipar, &
    ipivfl,ipivrb,isotri,iwrite,jhi,jlo, &
    ldafl,maxcofrb,maxelm,maxnew,maxnfl,maxnp, &
    maxny,maxopt,maxpar,maxparb,maxparf,maxsim, &
    nbcrb,ncofrb,nelem,neqnfl,nferb, &
    nlband,node,nodelm,np,npar,nparb,nparf,npe,nprof,nsenfl,ntay, &
    numnew,numopt,numsim,nx, &
    ny,par,parafl,pararb,pardif,parrb,parsav,parsen,partar, &
    phifl,phirb,rbase,rb,region,resfl,resflsav,resrb,reynld, &
    reytay,senfl, &
    senrb,splbmp,splflo,taubmp,tauflo,tecfil,tolnew,tolopt, &
    tolsim,value,wateb,watep,wateu,watev,wquad,xbl,xbr, &
    xc,xmax,xmin,xprof,xquad,xrange,xsiq, &
    ybl,ybr,yc,ymax,ymin,yquad,yrange)
!
!  Read the next command from the user
!
  command = '?'

  do

  if ( command(1:1) /= '#' .and. len_trim ( command ) /= 0 ) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Enter command:'
    end if

  end if

  read ( *, '(a)', iostat = ios ) command

  if ( ios /= 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ARBY4 - Fatal error!'
    write ( *, * ) '  I/O error for the input file.'
    stop
  end if

  if ( echo ) then
    write ( *, '(a)' ) trim ( command )
  end if
  write ( 17, '(a)' ) trim ( command )
!
!  Check for output-suppressing semicolon at end.
!
  if ( 0 < lenc ) then
    if ( command(lenc:lenc) == ';' ) then
      iprint = 0
      command(lenc:lenc) = ' '
      lenc = lenc-1
    else
      iprint = 1
    end if
  end if

15    continue

  if ( command(1:1) == '#') then
    cycle
  else if ( command == ' ') then
    cycle
  end if

  if ( 0 < iprint ) then
    write ( *, '(a)' ) ' '
  end if
!
!  COMPARE
!
  if ( s_eqi ( command,'compare')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - Compare:'
      write ( *, '(a)' ) '  Compare full solutions GFL and GFLSAV.'
    end if

    call fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar,nparf, &
      par,phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)

    call fxfl(area,eqn,gflsav,ifs,indx,nelem,neqnfl,node,np,npar, &
      nparf,par,phifl,region,resflsav,splflo,tauflo,xrange,yc,yrange)

    gfltmp(1:neqnfl) = gfl(1:neqnfl) - gflsav(1:neqnfl)
    resfltmp(1:neqnfl) = resfl(1:neqnfl) - resflsav(1:neqnfl)

    call nrmflo ( gfltmp, indx, neqnfl, np, resfltmp )
!
!  COSTFL
!
  else if ( s_eqi ( command,'costfl')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - Cost FL:'
      write ( *, '(a)' ) '  Evaluate the cost function for GFL.'
      write ( *, '(a)' ) ' '
    end if

    call getcst(cost,costb,costp,costu,costv,gfl,gfltar,indx,neqnfl,np, &
      nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)

    write ( *, '(a,g14.6)' ) 'Cost of GFL:',cost
!
!  COSTRB
!
  else if ( s_eqi ( command,'costrb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - Cost RB:'
      write ( *, '(a)' ) '  Evaluate the cost function for GRB.'
      write ( *, '(a)' ) ' '
    end if

    call grb2fl(gfl,gflrb,grb,maxnfl,ncofrb,neqnfl,rb)

    call getcst(cost,costb,costp,costu,costv,gfl,gfltar,indx,neqnfl,np, &
      nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)

    write ( *, '(a,g14.6)' ) 'Cost of GRB:',cost
!
!  DETFPFL
!
  else if ( s_eqi ( command,'detfpfl')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - DET FP FL'
      write ( *, '(a)' ) '  Compute determinant of full jacobian.'
    end if

    call fpfl(afl,area,eqn,gfl,indx,ldafl,maxelm,nelem,neqnfl, &
      nlband,node,np,npar,par,phifl)

    call dfacfl(afl,ldafl,neqnfl,nlband,nlband,ipivfl,info)

    call ddetfl(afl,detlog,detman,ipivfl,ldafl,neqnfl,nlband,nlband)

    write ( *, '(a,g14.6,a,g14.6)' ) 'Determinant = ', detman, &
      ' * 10 ** ', detlog
!
!  DETFPRB
!
  else if ( s_eqi ( command,'detfprb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - DET FP RB'
      write ( *, '(a)' ) '  Compute the reduced jacobian,'
      write ( *, '(a)' ) '  and its determinant.'
    end if

    call fprb(arb,area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb, &
      ncofrb,nelem,nferb,node,np,nx,ny,phirb,rb,reynld,xc,xrange,yc,yrange)

    call dfacrb(arb,maxcofrb,ncofrb,ipivrb,info)

    call ddetrb(arb,detlog,detman,ipivrb,maxcofrb,ncofrb)

    write ( *, '(a,g14.6,a,g14.6)' ) 'Determinant = ', detman, &
      ' * 10 ** ', detlog
!
!  DIFFPRB
!
  else if ( s_eqi ( command,'diffprb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - DIF FP RB'
      write ( *, '(a)' ) '  Estimate reduced jacobian by finite diff.'
    end if

    pararb(1:npar) = par(1:npar)
    grbarb(1:ncofrb) = grb(1:ncofrb)
    arb(1:maxcofrb,1:maxcofrb) = 0.0D+00

    call diffprb(arb,area,epsdif,grb,indx,maxcofrb,maxelm, &
      maxnfl,nbcrb,ncofrb,nelem,nferb,node,np,npar,nparf,nx, &
      ny,par,phirb,rb,resrb,reynld,tauflo,xc,xrange,yc,yrange)
!
!  DIFSENFL
!
  else if ( s_eqi ( command,'difsenfl')) then

    if ( ipar <= 0 .or. npar < ipar ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ARBY4 - Warning!'
      write ( *, '(a)' ) '  Cancelling the DIFSENFL command.'
      write ( *, '(a,i6)' ) '  IPAR = ',ipar
      write ( *, '(a)' ) '  but IPAR must be at least 0'
      write ( *, '(a,i6)' ) '  and no more than NPAR = ',npar
      cycle
    end if

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - DifSenFL:'
      write ( *, '(a)' ) '  Estimate full solution sensitivity '
      write ( *, '(a)' ) '  with respect to parameter ',ipar
      write ( *, '(a)' ) '  via finite differences.'
    end if

    ipar = npar
    dpar = drey

    parsen(1:npar) = par(1:npar)
    gflsen(1:neqnfl) = gfl(1:neqnfl)

    call difsenfl(afl,area,difcof,dpar,eqn,gfl,gflafl,ifs,ijac,indx,ipar, &
      ipivfl,iwrite,ldafl,maxcofrb,maxelm,maxnew,maxnfl,ncofrb,nelem,neqnfl, &
      nlband,node,np,npar,nparf,par,parafl,phifl,region,resfl,senfl,splflo, &
      tauflo,tolnew,xrange,yc,yrange)
!
!  DIFSENRB
!
  else if ( s_eqi ( command,'difsenrb')) then

    if ( ipar <= 0 .or. npar < ipar ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ARBY4 - Warning!'
      write ( *, '(a)' ) '  Cancelling the DIFSENRB command.'
      write ( *, '(a,i6)' ) '  IPAR = ', ipar
      write ( *, '(a)' ) '  but IPAR must be at least 0'
      write ( *, '(a,i6)' ) '  and no more than NPAR = ', npar
      cycle
    end if

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - Dif Sen RB:'
      write ( *, '(a)' ) '  Estimate reduced solution sensitivity '
      write ( *, '(a,i6)' ) '  with respect to parameter ', ipar
      write ( *, '(a)' ) '  via finite differences.'
    end if

    ipar = npar
    dpar = drey

    parsen(1:npar) = par(1:npar)
    grbsen(1:ncofrb) = grb(1:ncofrb)

    call difsenrb(arb,area,difcof,dpar,grb,grbarb,indx,ipar,ipivrb,iwrite, &
      maxcofrb,maxelm,maxnew,maxnfl,nbcrb,ncofrb,nelem,nferb,node,np,npar, &
      nparf,nx,ny,par,pararb,phirb,rb,resrb,senrb,tauflo,tolnew,xc,xrange, &
      yc,yrange)
!
!  DisPlot
!
  else if ( s_eqi ( command,'displot')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - DisPlot: Write data to DISPLAY plot file.'
    end if

    call intprs(gfl,indx,nelem,neqnfl,node,np,p)

    u(1:np) = gfl(indx(1,1:np))
    v(1:np) = gfl(indx(2,1:np))

    call wrdis(disfil,eqn,indx,isotri,maxcofrb,maxnfl,ncofrb,nelem,neqnfl, &
      node,np,npar,npe,nprof,nx,ny,p,par,rb,u,v,xc,xprof,yc)
!
!  DREY =
!
  else if ( s_eqi ( command(1:4),'drey')) then

    if ( s_eqi ( command(1:5),'drey=')) then
      read(command(6:),*)drey
    else
      write ( *, '(a)' ) 'Enter value for DREY:'
      read(*,*)drey
      write(17,*)drey
    end if

    if ( 0 < iprint ) then
      write ( *, '(a,g14.6)' ) 'DREY set to ',drey
    end if
!
!  EPSDIF =
!
  else if ( s_eqi ( command(1:6),'epsdif')) then

    if ( s_eqi ( command(1:7),'epsdif=')) then
      read(command(8:),*)epsdif
    else
      write ( *, '(a)' ) 'Enter value for EPSDIF:'
      read ( *, * ) epsdif
      write ( 17, '(g14.6)' ) epsdif
    end if

    if ( 0 < iprint ) then
      write ( *, '(a,g14.6)' ) 'EPSDIF set to ', epsdif
    end if
!
!  ECHO
!
  else if ( s_eqi ( command,'echo')) then

    echo = .not.echo
    if ( echo) then
      write(*,'(a)')command
      if ( 0 < iprint ) then
        write ( *, '(a)' ) 'User commands will be echoed.'
      end if
    else
      if ( 0 < iprint ) then
        write ( *, '(a)' ) 'User commands will not be echoed.'
      end if
    end if
!
!  EXPAND GRB
!
  else if ( s_eqi ( command,'expand grb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - Expand GRB:'
      write ( *, '(a)' ) '  Compute GFL = RB * GRB'
    end if

   call grb2fl(gfl,gflrb,grb,maxnfl,ncofrb,neqnfl,rb)
!
!  DISFIL =
!
  else if ( s_eqi ( command(1:6),'disfil')) then

    if ( s_eqi ( command(1:7),'disfil=')) then
      disfil = command(8:)
    else
      write ( *, '(a)' ) 'Enter the DISPLAY output file name:'
      read(*,'(a)')disfil
      write(17,'(a)')disfil
    end if

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'DISPLAY output file name set to ' // trim ( disfil )
    end if
!
!  FILTEC =
!
  else if ( s_eqi ( command(1:6),'tecfil')) then

    if ( s_eqi ( command(1:7),'tecfil=')) then
      tecfil = command(8:)
    else
      write ( *, '(a)' ) 'Enter the TECPLOT output file name:'
      read(*,'(a)')tecfil
      write(17,'(a)')tecfil
    end if

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'TECPLOT output file name set to ' // trim ( tecfil )
    end if
!
!  FPFL
!
  else if ( s_eqi ( command,'fpfl')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FPFL - Evaluate full jacobian.'
    end if

    call fpfl(afl,area,eqn,gfl,indx,ldafl,maxelm,nelem,neqnfl, &
      nlband,node,np,npar,par,phifl)
!
!  FPIRB
!
  else if ( s_eqi ( command,'fpirb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FPIRB:'
      write ( *, '(a)' ) '  Evaluate FP indirectly at reduced solution GRB,'
      write ( *, '(a)' ) '  by expanding GRB to GFL, evaluating FP(GFL),'
      write ( *, '(a)' ) '  and then reducing.'
    end if

    call fpirb(afl,arb,area,eqn,gflrb,grb,indx,ldafl,maxcofrb, &
      maxelm,maxnfl,nbcrb,ncofrb,nelem,neqnfl,nferb,nlband,node, &
      np,npar,nx,ny,par,phifl,rb,xc,xrange,yc,yrange)
!
!  FPRB
!
  else if ( s_eqi ( command,'fprb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FPRB - Evaluate reduced jacobian.'
    end if

    pararb(1:npar) = par(1:npar)
    grbarb(1:ncofrb) = grb(1:ncofrb)

    call fprb(arb,area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb,nelem, &
      nferb,node,np,nx,ny,phirb,rb,reynld,xc,xrange,yc,yrange)
!
!  FPRB = 0
!
  else if ( s_eqi ( command,'fprb=0')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FPRB = 0 - Zero out the reduced jacobian.'
    end if

    arb(1:maxcofrb,1:maxcofrb) = 0.0D+00
!
!  FXFL
!
  else if ( s_eqi ( command,'fxfl')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FX FL:'
      write ( *, '(a)' ) '  Evaluate FXFL at full solution GFL.'
    end if

    call fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar,nparf,par, &
      phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)
!
!  FXIRB
!
  else if ( s_eqi ( command,'fxirb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FX IRB:'
      write ( *, '(a)' ) '  Evaluate FX indirectly at reduced solution GRB,'
      write ( *, '(a)' ) '  by expanding GRB to GFL, evaluating F(GFL),'
      write ( *, '(a)' ) '  and then reducing.'
    end if

    call fxirb(area,eqn,gflrb,grb,ifs,indx,maxcofrb,maxnfl,nbcrb,ncofrb, &
      nelem,neqnfl,nferb,node,np,npar,nparf,nx,ny,par,phifl,rb,region,resrb, &
      splflo,tauflo,xc,xrange,yc,yrange)
!
!  FXRB
!
  else if ( s_eqi ( command,'fxrb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FX RB:'
      write ( *, '(a)' ) '  Evaluate FX directly at reduced solution GRB.'
    end if

    reynld = par(npar)

    call fxrb(area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb,nelem,nferb, &
      node,np,npar,nparf,nx,ny,par,phirb,rb,resrb,reynld,tauflo,xc,xrange, &
      yc,yrange)
!
!  FXRB = 0
!
  else if ( s_eqi ( command,'fxrb=0')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - FXRB = 0'
    end if

    resrb(1:ncofrb) = 0.0D+00
!
!  GETGSEN
!
  else if ( s_eqi ( command,'getgsen')) then

    call getgsen(grb,gsen,icolrb,maxcofrb,nbcrb,ncofrb,nsenfl,rbase)
!
!  GETRB
!
  else if ( s_eqi ( command,'getrb')) then

    if ( 0 < iprint ) then
      write ( *, '(a)' ) 'ARBY4 - Get RB:'
      write ( *, '(a)' ) '  Compute reduced basis RB at current GFL.'
    end if

    if ( ncofrb < 0) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ARBY4 - Warning!'
      write ( *, '(a)' ) '  The GETRB command is being cancelled,'
      write ( *, '(a,i6)' ) '  since NCOFRB is ', ncofrb
      write ( *, '(a)' ) '  Please use the "NCOFRB = " command first!'
    end if

    if ( dvneq ( npar, par, parsen ) ) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Error!'
      write ( *, * ) '  Please compute the sensitivities first,'
      write ( *, * ) '  using the GETSENFL command!'
      cycle
    end if
!
!  Save the affine displacement vector GFLRB.
!
    parrb(1:npar) = par(1:npar)
    gflrb(1:neqnfl) = gfl(1:neqnfl)
!
!  Get boundary condition reduced basis vectors.
!
    call getbcrb(gflrb,maxcofrb,maxnfl,nbcrb,neqnfl,rb)
!
!  Get finite element reduced basis vectors.
!
    call getferb(icolrb,maxcofrb,maxnfl,nbcrb,ncofrb,neqnfl, &
      nferb,nsenfl,rb,rbase,senfl,senrb)

    grb(1:ncofrb) = 0.0D+00
!
!  Compute PHIRB, the reduced basis functions at quadrature points.
!
    if ( 0 < iprint ) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Note:'
      write ( *, * ) '  Automatically evaluating PHIRB, the reduced'
      write ( *, * ) '  basis functions at quadrature points.'
    end if

    call setprb(eqn,indx,maxcofrb,maxelm,maxnfl,nelem,neqnfl, &
      ncofrb,node,np,phifl,phirb,rb)
!
!  GETSENFL
!
  else if ( s_eqi ( command,'getsenfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Get Sen FL:'
      write ( *, * ) '  Compute full solution sensitivities with '
      write ( *, * ) '  respect to the parameter REYNLD, up to order'
      write ( *, * ) '  NSENFL = ',nsenfl
    end if

    if ( nsenfl < 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  The GETSENFL command is being cancelled,'
      write ( *, * ) '  since NSENFL is ',nsenfl
      write ( *, * ) '  Please use the "NSENFL = " command first!'
    else

      parsen(1:npar) = par(1:npar)
      gflsen(1:neqnfl) = gfl(1:neqnfl)

      call getsenfl(afl,area,eqn,gfl,indx,ipivfl,ldafl,maxcofrb,maxnfl, &
        nelem,neqnfl,nlband,node,np,npar,nsenfl,par,phifl,resfl,senfl)

    end if
!
!  GETSENRB
!
  else if ( s_eqi ( command,'getsenrb')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Get Sen RB'
      write ( *, * ) '  Compute reduced basis sensitivities.'
    end if

    call getsenrb(maxcofrb,maxnfl,ncofrb,neqnfl,nsenfl,rb,senfl,senrb)
!
!  GFL = 0, GFLSAV, GFLTAY, TAYLOR
!
  else if ( s_eqi ( command(1:4),'gfl=')) then

    if ( command(5:5) == '0') then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFL = 0'
        write ( *, * ) '  Set full solution estimate GFL to zero.'
      end if

      gfl(1:neqnfl) = 0.0D+00

    else if ( s_eqi ( command(5:10),'gflsav')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFL = GFLSAV'
        write ( *, * ) '  Set full solution estimate GFL to GFLSAV.'
      end if

      gfl(1:neqnfl) = gflsav(1:neqnfl)

    else if ( s_eqi ( command(5:10),'gfltay')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFL = GFLTAY'
        write ( *, * ) '  Set full solution estimate GFL to GFLTAY.'
      end if

      gfl(1:neqnfl) = gfltay(1:neqnfl)

    else if ( s_eqi ( command(5:10),'taylor')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFL = TAYLOR'
        write ( *, * ) '  Set full solution GFL to Taylor prediction,'
        write ( *, * ) '  based at REYTAY = ',reytay
        write ( *, * ) '  using NTAY = ',ntay,' terms.'
      end if

      gfl(1:neqnfl) = gfltay(1:neqnfl)

      do i = 1, neqnfl
        do j = 1, ntay
          call fact(j,factj)
          temp = ((reynld-reytay)**j)/factj
          gfl(i) = gfl(i)+temp*senfl(i,j)
        end do
      end do

    end if
!
!  GFLSAV = GFL
!
  else if ( s_eqi ( command,'gflsav=gfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - GFLSAV = GFL'
      write ( *, * ) '  Save current full solution.'
    end if

    parsav(1:npar) = par(1:npar)
    gflsav(1:neqnfl) = gfl(1:neqnfl)
    resflsav(1:neqnfl) = resfl(1:neqnfl)
!
!  GFLTAY = 0, GFL, GFLSAV
!
  else if ( s_eqi ( command(1:7),'gfltay=')) then

    if ( command(8:8) == '0') then
      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTAY = 0'
        write ( *, * ) '  Set Taylor base full solution to zero.'
      end if

      gfltay(1:neqnfl) = 0.0D+00

    else if ( s_eqi ( command(8:13),'gflsav')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTAY = GFLSAV'
        write ( *, * ) '  Set Taylor base full solution to GFLSAV.'
      end if

      gfltay(1:neqnfl) = gflsav(1:neqnfl)

    else if ( s_eqi ( command(8:10),'gfl')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTAY = GFL'
        write ( *, * ) '  Set Taylor base full solution to GFL.'
      end if

      gfltay(1:neqnfl) = gfl(1:neqnfl)

    end if
!
!  GFLTMP = 0/GFL/GFL-GFLSAV/GFL-GFLTAR/GFLSAV/GFLSAV-GFLTAY
!
  else if ( s_eqi ( command(1:7),'gfltmp=')) then

    if ( command(8:8) == '0') then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTMP = 0'
      end if

      gfltmp(1:neqnfl) = 0.0D+00

    else if ( s_eqi ( command(8:),'gfl')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTMP = GFL'
      end if

      gfltmp(1:neqnfl) = gfl(1:neqnfl)

    else if ( s_eqi ( command(8:),'gfl-gflsav')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTMP = GFL-GFLSAV'
      end if

      gfltmp(1:neqnfl) = gfl(1:neqnfl) - gflsav(1:neqnfl)

    else if ( s_eqi ( command(8:),'gfl-gfltar')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTMP = GFL-GFLTAR'
      end if

      gfltmp(1:neqnfl) = gfl(1:neqnfl) - gfltar(1:neqnfl)

    else if ( s_eqi ( command(8:),'gflsav')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTMP = GFLSAV'
      end if

      gfltmp(1:neqnfl) = gflsav(1:neqnfl)

    else if ( s_eqi ( command(8:),'gflsav-gfltay')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GFLTMP = GFLSAV-GFLTAY'
      end if

      gfltmp(1:neqnfl) = gflsav(1:neqnfl) - gfltay(1:neqnfl)

    else

      write ( *, * ) 'ARBY4 - Error'
      write ( *, * ) '  Did not understand your command!'
      cycle

    end if
!
!  GRB(*) = *
!
  else if ( s_eqi ( command(1:4),'grb(')) then

    call chrcti(command(5:),ival,ierror,lchar)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  ChrCTI returned nonzero error flag!'
      cycle
    end if

    if ( ival < 1 .or. ncofrb < ival ) then
      write ( *, * ) ' '
      write ( *, * ) 'INPUT - Warning!'
      write ( *, * ) '  Index IVAL of GRB is out of bounds!'
      write ( *, * ) '  IVAL = ',ival
      cycle
    end if

    call chrctd(command(5+lchar+2:),value,ierror,lchar)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  ChrCTD returned nonzero error flag!'
      cycle
    end if

    grb(ival) = value

    if ( 0 < iprint ) then
      write ( *, * ) 'GRB(',ival,') set to ',grb(ival)
    end if
!
!  GRB = 0, GRBSAV, GRBTAY, TAYLOR
!
  else if ( s_eqi ( command(1:4),'grb=')) then

    if ( s_eqi ( command(1:5),'grb=0')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRB = 0'
        write ( *, * ) '  Set reduced solution to 0.'
      end if

      grb(1:ncofrb) = 0.0D+00

    else if ( s_eqi ( command,'grb=grbsav')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRB = GRBSAV'
        write ( *, * ) '  Set reduced solution to saved value.'
      end if

      grb(1:ncofrb) = grbsav(1:ncofrb)

    else if ( s_eqi ( command,'grb=grbtay')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRB = GRBTAY'
        write ( *, * ) '  Set reduced solution to Taylor base.'
      end if

      grb(1:ncofrb) = grbtay(1:ncofrb)

    else if ( s_eqi ( command,'grb=taylor')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRB = TAYLOR'
        write ( *, * ) '  Set reduced solution to Taylor prediction,'
        write ( *, * ) '  based at REYTAY, GRBTAY plus sensitivities.'
      end if

      do i = 1, ncofrb
        grb(i) = 0.0D+00
        do j = 1, ntay
          jtay = j-1
          call fact(jtay,factj)
          temp = ((reynld-reytay)**jtay)/factj
          grb(i) = grb(i)+temp*senrb(i,j)
        end do
      end do

    else if ( s_eqi ( command(1:5),'grb=(')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRB = (v0,v1,...,vNeqnRB)'
      end if

      klo = 6
      do i = 1, ncofrb
        call chrctd(command(klo:),temp,ierror,lchar)
        if ( ierror /= 0) then
          write ( *, * ) 'ARBY4 - Warning!'
          write ( *, * ) '  There was an error reading your data.'
          grb(i) = 0.0D+00
        else
          grb(i) = temp
        end if
        klo = klo+lchar
      end do

    end if
!
!  GRBSAV = 0, GRB
!
  else if ( s_eqi ( command(1:7),'grbsav=')) then

    if ( command(8:8) == '0') then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRBSAV = 0'
      end if

      grbsav(1:ncofrb) = 0.0D+00

    else if ( s_eqi ( command,'grbsav=grb')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRBSAV = GRB'
      end if

      grbsav(1:ncofrb) = grb(1:ncofrb)

    end if
!
!  GRBTAY = 0, GRB, GRBSAV
!
  else if ( s_eqi ( command(1:7),'grbtay=')) then

    if ( command(8:8) == '0') then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRBTAY = 0'
        write ( *, * ) '  Set Taylor base reduced solution to zero.'
      end if

      grbtay(1:ncofrb) = 0.0D+00

    else if ( s_eqi ( command(8:13),'grbsav')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRBTAY = GRBSAV'
        write ( *, * ) '  Set Taylor base reduced solution to GRBSAV.'
      end if

      grbtay(1:ncofrb) = grbsav(1:ncofrb)

    else if ( s_eqi ( command(8:10),'grb')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - GRBTAY = GRB'
        write ( *, * ) '  Set Taylor base reduced solution to GRB.'
      end if

      grbtay(1:ncofrb) = grb(1:ncofrb)

    end if
!
!  GRIDX =
!
  else if ( s_eqi ( command(1:5),'gridx')) then

    if ( s_eqi ( command(1:6),'gridx=')) then
      gridx = command(7:)
    else
      write ( *, * ) 'Enter GRIDX option: UNIFORM, COS, SINSQ:'
      read(*,'(a)')gridx
      write(17,'(a)')gridx
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'The GRIDX option set to '//gridx
      write ( *, * ) 'Remember to use the SETGEO command'
      write ( *, * ) 'before trying to solve your system!'
    end if
!
!  GRIDY =
!
  else if ( s_eqi ( command(1:5),'gridy')) then

    if ( s_eqi ( command(1:6),'gridy=')) then
      gridy = command(7:)
    else
      write ( *, * ) 'Enter GRIDY option: UNIFORM, COS, SINSQ:'
      read(*,'(a)')gridy
      write(17,'(a)')gridy
    end if

    if ( 0 < iprint ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'The GRIDY option set to ' // trim ( gridy )
      write ( *, '(a)' ) 'Remember to use the SETGEO command'
      write ( *, '(a)' ) 'before trying to solve your system!'
    end if
!
!  HELLO
!
  else if ( s_eqi ( command,'hello')) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'ARBY4 - Hello'
    write ( *, '(a)' ) ' '

    call hello ( maxnx, maxny )
!
!  HELP
!
  else if ( s_eqi ( command,'help')) then

    call help
!
!  IBS =
!
  else if ( s_eqi ( command(1:3),'ibs')) then

    if ( s_eqi ( command(1:4),'ibs=')) then
      read(command(5:),*)ibs
    else
      write ( *, * ) 'Enter value for IBS:'
      read(*,*)ibs
      write(17,*)ibs
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'IBS set to ',ibs
    end if
!
!  IBUMP =
!
  else if ( s_eqi ( command(1:5),'ibump')) then

    if ( s_eqi ( command(1:6),'ibump=')) then
      read(command(7:),*)ibump
    else
      write ( *, * ) 'Enter value for IBUMP:'
      read(*,*)ibump
      write(17,*)ibump
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'IBUMP set to ',ibump
    end if
!
!  IFS =
!
  else if ( s_eqi ( command(1:3),'ifs')) then

    if ( s_eqi ( command(1:4),'ifs=')) then
      read(command(5:),*)ifs
    else
      write ( *, * ) 'Enter value for IFS:'
      read(*,*)ifs
      write(17,*)ifs
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'IFS set to ',ifs
    end if
!
!  IHI =
!
  else if ( s_eqi ( command(1:3),'ihi')) then

    if ( s_eqi ( command(1:4),'ihi=')) then
      if ( s_eqi ( command,'ihi=ncofrb')) then
        ihi = ncofrb
      else if ( s_eqi ( command,'ihi=neqnfl')) then
        ihi = neqnfl
      else if ( s_eqi ( command,'ihi=np')) then
        ihi = np
      else
        read(command(5:),*)ihi
      end if
    else
      write ( *, * ) 'Enter value for IHI:'
      read(*,*)ihi
      write(17,*)ihi
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'IHI set to ',ihi
    end if
!
!  IJAC =
!
  else if ( s_eqi ( command(1:4),'ijac')) then

    if ( s_eqi ( command(1:5),'ijac=')) then
      read(command(6:),*)ijac
    else
      write ( *, * ) 'Enter value for IJAC:'
      read(*,*)ijac
      write(17,*)ijac
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'IJAC set to ',ijac
    end if
!
!  ILO =
!
  else if ( s_eqi ( command(1:3),'ilo')) then

    if ( s_eqi ( command(1:4),'ilo=')) then
      read(command(5:),*)ilo
    else
      write ( *, * ) 'Enter value for ILO:'
      read(*,*)ilo
      write(17,*)ilo
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'ILO set to ',ilo
    end if
!
!  INIT
!
  else if ( s_eqi ( command,'init')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Init'
      write ( *, * ) '  Initialize all data to zero.'
      write ( *, * ) ' '
    end if

    call init(afl,arb,area,command,cost,costb,costp,costu,costv,difcof, &
      disfil,drey,epsdif,eqn,etaq,gfl,gflafl,gflrb,gflsav,gflsen,gfltar, &
      gfltay,grb,grbarb,grbsav,grbsen,grbtay,gridx,gridy,hx,hy,ibs,ibump, &
      icolrb,ierror,ifs,ihi,ijac,ilo,indx,iopt,ipar,ipivfl,ipivrb,isotri, &
      iwrite,jhi,jlo,ldafl,maxcofrb,maxelm,maxnew,maxnfl,maxnp,maxny, &
      maxopt,maxpar,maxparb,maxparf,maxsim,nbcrb,ncofrb,nelem,neqnfl, &
      nferb,nlband,node,nodelm,np,npar,nparb,nparf,npe,nprof,nsenfl,ntay, &
      numnew,numopt,numsim,nx,ny,par,parafl,pararb,pardif,parrb,parsav, &
      parsen,partar,phifl,phirb,rbase,rb,region,resfl,resflsav,resrb, &
      reynld,reytay,senfl,senrb,splbmp,splflo,taubmp,tauflo,tecfil,tolnew, &
      tolopt,tolsim,value,wateb,watep,wateu,watev,wquad,xbl,xbr,xc, &
      xmax,xmin,xprof,xquad,xrange,xsiq,ybl,ybr,yc,ymax,ymin,yquad,yrange)
!
!  IOPT(*) =
!
    else if ( s_eqi ( command(1:5),'iopt(')) then

      call chrcti(command(6:),ival1,ierror,lchar)

      if ( ierror /= 0) then
        write ( *, * ) ' '
        write ( *, * ) 'INPUT - Fatal error!'
        write ( *, * ) '  ChrCTI returned nonzero error flag!'
        stop
      end if

      if ( ival1 < 1 .or. maxpar < ival1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'INPUT - Fatal error!'
        write ( *, * ) '  Index IVAL1 of IOPT is out of bounds!'
        write ( *, * ) '  IVAL1 = ', ival1
        stop
      end if

      call chrcti(command(6+lchar+2:),ival2,ierror,lchar)

      iopt(ival1) = ival2

      if ( 0 < iprint ) then
        write ( *, * ) 'IOPT(',ival1,' ) set to ',ival2
      end if
!
!  IWRITE =
!
  else if ( s_eqi ( command(1:6),'iwrite')) then

    if ( s_eqi ( command(1:7),'iwrite=')) then
      read(command(8:),*)iwrite
    else
      write ( *, * ) 'Enter value for IWRITE:'
      read(*,*)iwrite
      write(17,*)iwrite
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'IWRITE set to ',iwrite
    end if
!
!  JHI =
!
  else if ( s_eqi ( command(1:3),'jhi')) then

    if ( s_eqi ( command(1:4),'jhi=')) then
      if ( s_eqi ( command,'jhi=ncofrb')) then
        jhi = ncofrb
      else if ( s_eqi ( command,'jhi=neqnfl')) then
        jhi = neqnfl
      else if ( s_eqi ( command,'jhi=nsenfl')) then
        jhi = nsenfl
      else
        read(command(5:),*)jhi
      end if
    else
      write ( *, * ) 'Enter value for JHI:'
      read(*,*)jhi
      write(17,*)jhi
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'JHI set to ',jhi
    end if
!
!  JLO =
!
  else if ( s_eqi ( command(1:3),'jlo')) then

    if ( s_eqi ( command(1:4),'jlo=')) then
      read(command(5:),*)jlo
    else
      write ( *, * ) 'Enter value for JLO:'
      read(*,*)jlo
      write(17,*)jlo
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'JLO set to ',jlo
    end if
!
!  L2NORM GFL/GFLSAV/GFLTAR/GFLTAY/GFLTMP
!
  else if ( s_eqi ( command(1:6),'l2norm')) then

    if ( s_eqi ( command(8:),'gfl')) then
      call l2norm(gfl,gflnrm,indx,nelem,neqnfl,node,np,xc,yc)
      write ( *, * ) 'ARBY4 - L2Norm of GFL = ',gflnrm
    else if ( s_eqi ( command(8:),'gflsav')) then
      call l2norm(gflsav,gflnrm,indx,nelem,neqnfl,node,np,xc,yc)
      write ( *, * ) 'ARBY4 - L2Norm of GFLSAV = ',gflnrm
    else if ( s_eqi ( command(8:),'gfltar')) then
      call l2norm(gfltar,gflnrm,indx,nelem,neqnfl,node,np,xc,yc)
      write ( *, * ) 'ARBY4 - L2Norm of GFLTAR = ',gflnrm
    else if ( s_eqi ( command(8:),'gfltay')) then
      call l2norm(gfltay,gflnrm,indx,nelem,neqnfl,node,np,xc,yc)
      write ( *, * ) 'ARBY4 - L2Norm of GFLTAY = ',gflnrm
    else if ( s_eqi ( command(8:),'gfltmp')) then
      call l2norm(gfltmp,gflnrm,indx,nelem,neqnfl,node,np,xc,yc)
      write ( *, * ) 'ARBY4 - L2Norm of GFLTMP = ',gflnrm
    else
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Error!'
      write ( *, * ) '  Legal choices were GFL/GFLSAV/GFLTAY/GFLTMP.'
      write ( *, * ) '  Your choice was '//command(8:)
      cycle
    end if
!
!  MAXNEW =
!
  else if ( s_eqi ( command(1:6),'maxnew')) then

    if ( s_eqi ( command(1:7),'maxnew=')) then
      read(command(8:),*)maxnew
    else
      write ( *, * ) 'Enter value for MAXNEW:'
      read(*,*)maxnew
      write(17,*)maxnew
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'MAXNEW set to ',maxnew
    end if
!
!  MAXOPT =
!
  else if ( s_eqi ( command(1:6),'maxopt')) then

    if ( s_eqi ( command(1:7),'maxopt=')) then
      read(command(8:),*)maxopt
    else
      write ( *, * ) 'Enter value for MAXOPT:'
      read(*,*)maxopt
      write(17,*)maxopt
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'MAXOPT set to ',maxopt
    end if
!
!  MAXSIM =
!
  else if ( s_eqi ( command(1:6),'maxsim')) then

    if ( s_eqi ( command(1:7),'maxsim=')) then
      read(command(8:),*)maxsim
    else
      write ( *, * ) 'Enter value for MAXSIM:'
      read(*,*)maxsim
      write(17,*)maxsim
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'MAXSIM set to ',maxsim
    end if
!
!  NBCRB =
!
  else if ( s_eqi ( command(1:5),'nbcrb')) then

    if ( s_eqi ( command(1:6),'nbcrb=')) then
      read(command(7:),*)nbcrb
    else
      write ( *, * ) 'Enter value for NBCRB:'
      read(*,*)nbcrb
      write(17,*)nbcrb
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'NBCRB set to ',nbcrb
    end if
!
!  NEWTFL
!  Apply Newton's method to full solution estimate.
!
  else if ( s_eqi ( command,'newtfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - NewtFL'
      write ( *, * ) '  Apply Newton to full solution estimate GFL.'
    end if

    call newtfl(afl,area,eqn,gfl,gflafl,ierror,ifs,ijac,indx,ipivfl,iwrite, &
      ldafl,maxelm,maxnew,nelem,neqnfl,nlband,node,np,npar,nparf,numnew,par, &
      parafl,phifl,region,resfl,rmax,splflo,tauflo,tolnew,xrange,yc,yrange)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Fatal error!'
      write ( *, * ) '  NEWTFL failed!'
      write ( *, * ) '  The parameters at which failure occurred:'
      write ( *, * ) ' '
      call prpar(iopt,npar,nparb,nparf,par)
      ierror = 1
      cycle
    else
      if ( iwrite <= 1) then
        write ( *, * ) '  Newton step ',numsim,' residual norm  = ',rmax
      end if
    end if
!
!  NEWTRB
!  Apply Newton's method to reduced solution estimate.
!
  else if ( s_eqi ( command,'newtrb')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - NewtRB'
      write ( *, * ) '  Apply Newton to reduced solution estimate GRB.'
    end if

    if ( ncofrb <= 0 ) then
      cycle
    end if

    call newtrb(arb,area,grb,grbarb,ierror,indx,ipivrb, &
      iwrite,maxcofrb,maxelm,maxnew,maxnfl,nbcrb,ncofrb,nelem, &
      nferb,node,np,npar,nparf,nx,ny,par,pararb,phirb, &
      rb,resrb,rmax,tauflo,tolnew,xc,xrange,yc,yrange)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Fatal error!'
      write ( *, * ) '  NEWTRB failed!'
      write ( *, * ) '  The parameters at which failure occurred:'
      call prpar(iopt,npar,nparb,nparf,par)
      ierror = 1
      cycle
    else
      if ( iwrite <= 1) then
        write ( *, * ) '  Final Newton residual was MxNorm(FXRB) = ',rmax
      end if
    end if
!
!  NPARB =
!
  else if ( s_eqi ( command(1:5),'nparb')) then

    if ( s_eqi ( command(1:6),'nparb=')) then
      read(command(7:),*)nparb
    else
      write ( *, * ) 'Enter value for NPARB:'
      read(*,*)nparb
      write(17,*)nparb
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'NPARB set to ',nparb
    end if
!
!  NPARF =
!
  else if ( s_eqi ( command(1:5),'nparf')) then

    if ( s_eqi ( command(1:6),'nparf=')) then
      read(command(7:),*)nparf
    else
      write ( *, * ) 'Enter value for NPARF:'
      read(*,*)nparf
      write(17,*)nparf
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'NPARF set to ',nparf
    end if
!
!  NSENFL = #
!
  else if ( s_eqi ( command(1:6),'nsenfl')) then

    if ( s_eqi ( command(1:7),'nsenfl=')) then
      read(command(8:),*)itemp
    else
      write ( *, * ) 'Enter value for NSENFL:'
      read(*,*)itemp
      write(17,*)itemp
    end if

    if ( itemp < 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  NSENFL must be at least 0!'
      write ( *, * ) '  but your value was ',itemp
    else
      nsenfl = itemp
      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - NSENFL set to ',nsenfl
      end if
    end if
!
!  NTAY = #
!  NTAY = NCOFRB
!
  else if ( s_eqi ( command(1:4),'ntay')) then

    if ( s_eqi ( command(1:5),'ntay=')) then
      if ( s_eqi ( command(6:11),'ncofrb')) then
        itemp = ncofrb
      else
        read(command(6:),*)itemp
      end if
    else
      write ( *, * ) 'Enter value for NTAY:'
      read(*,*)itemp
      write(17,*)itemp
    end if

    if ( itemp < 0 .or. ncofrb < itemp ) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  NTAY must be between 0 and NCOFRB =  ',ncofrb
      write ( *, * ) '  but your value was ',itemp
    else
      ntay = itemp
      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - NTAY set to ',ntay
      end if
    end if
!
!  NX =
!
  else if ( s_eqi ( command(1:2),'nx')) then

    if ( s_eqi ( command(1:3),'nx=')) then
      read(command(4:),*)itemp
    else
      write ( *, * ) 'Enter value for NX:'
      read(*,*)itemp
      write(17,*)itemp
    end if

    if ( itemp < 2 ) then
      write ( *, * ) 'ARBY4 - Unacceptable input.'
      write ( *, * ) '  NX must be at least 2.'
      write ( *, * ) '  Your value was ',itemp
    else if ( maxnx < itemp ) then
      write ( *, * ) 'ARBY4 - Unacceptable input.'
      write ( *, * ) '  NX must be no more than MAXNX = ',maxnx
      write ( *, * ) '  Your value was ',itemp
    else
      nx = itemp
      if ( 0 < iprint ) then
        write ( *, * ) 'NX set to ',nx
        write ( *, * ) 'Remember to use the SETLOG and SETGEO commands'
        write ( *, * ) 'before trying to solve your systems!'
      end if
    end if
!
!  NY =
!
  else if ( s_eqi ( command(1:2),'ny')) then

    if ( s_eqi ( command(1:3),'ny=')) then
      read(command(4:),*)itemp
    else
      write ( *, * ) 'Enter value for NY:'
      read(*,*)itemp
      write(17,*)itemp
    end if

    if ( itemp < 2) then
      write ( *, * ) 'ARBY4 - Unacceptable input.'
      write ( *, * ) '  NY must be at least 2.'
      write ( *, * ) '  Your value was ',itemp
    else if ( maxny < itemp ) then
      write ( *, * ) 'ARBY4 - Unacceptable input.'
      write ( *, * ) '  NY must be no more than MAXNY = ',maxny
      write ( *, * ) '  Your value was ',itemp
    else
      ny = itemp
      if ( 0 < iprint ) then
        write ( *, * ) 'NY set to ',ny
        write ( *, * ) 'Remember to use the SETLOG and SETGEO commands'
        write ( *, * ) 'before trying to solve your system!'
      end if
    end if
!
!  OPTFL
!
  else if ( s_eqi ( command,'optfl')) then

    write ( *, * ) 'ARBY4 - OPT FL:'
    write ( *, * ) '  Optimize the full system.'
    write ( *, * ) '  NO WAY, JOSE!  NOT READY YET!'
!
!  OPTDIFFL
!
  else if ( s_eqi ( command,'optdiffl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - OptDifFl:'
      write ( *, * ) '  Optimize the cost of the full system;'
      write ( *, * ) '  The optimization code will approximate cost'
      write ( *, * ) '  gradients by finite differences.'
      write ( *, * ) '  Initial estimate is (PAR,GFL,COST).'
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Note!'
      write ( *, * ) '  You must already have issued the TARGET command!'
    end if

    call optdiffl(afl,area,cost,dopt,eqn,etaq,gfl,gflafl, &
         gflopt,gfltar,gridx,gridy,ibs,ierror,ifs,ijac,indx,iopt, &
         ipivfl,isotri,ivopt,iwrite,ldafl,liv,lv,maxelm,maxnew,maxnfl, &
         maxnp,maxny,maxopt,maxpar,maxparb,maxparf,maxsim,nelem, &
         neqnfl,nlband,node,nodelm,np,npar,nparb,nparf,nprof,numdif, &
         numopt,nx,ny,par,parafl,paropt,phifl,region,resfl,splbmp, &
         splflo,taubmp,tauflo,tolnew,tolopt,tolsim,vopt,wateb,watep, &
         wateu,watev,wquad,xbl,xbr,xc,xopt,xquad,xrange,xsiq,ybl, &
         ybr,yc,yquad,yrange)

    if ( ierror == 0) then

      par(1:npar) = paropt(1:npar)
      gfl(1:neqnfl) = gflopt(1:neqnfl)

      write ( *, * ) ' '
      write ( *, * ) 'Optimizing parameters:'
      call prpar(iopt,npar,nparb,nparf,par)
      write ( *, * ) ' '
      write ( *, * ) 'Optimal cost = ',cost
      write ( *, * ) ' '
      write ( *, * ) 'Number of standard full solutions:',numopt
      write ( *, * ) 'Number of auxilliary solutions:   ',numdif

    else

      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  The optimization was unsuccessful.'
    end if
!
!  OPTDIFRB
!
  else if ( s_eqi ( command,'optdifrb')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - OptDifRB:'
      write ( *, * ) '  Optimize the cost of the reduced system;'
      write ( *, * ) '  The optimization code will approximate cost'
      write ( *, * ) '  gradients by finite differences.'
      write ( *, * ) '  Initial estimate is (PAR,GRB,COST).'
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Note!'
      write ( *, * ) '  You must already have issued the commands:'
      write ( *, * ) '  SETLOG, SETGEO, NEWTFL, GETSEN, GETRB, TARGET!'
    end if

    call optdifrb(arb,area,cost,dopt,gflrb,gfltar,gfltmp, &
         grb,grbarb,grbopt,ierror,indx,iopt,ipivrb, &
         ivopt,iwrite,liv,lv,maxcofrb,maxelm,maxnew,maxnfl,maxnp, &
         maxny,maxopt,maxpar,maxparb,maxsim,nbcrb,ncofrb,nelem,neqnfl, &
         nferb,node,np,npar,nparb,nparf,nprof,numdif,numopt,nx,ny,par, &
         pararb,paropt,phirb,rb,resrb,splbmp,tauflo,taubmp,tolnew, &
         tolopt,tolsim,vopt,wateb,watep,wateu,watev,xbl,xbr,xc,xopt, &
         xrange,ybl,ybr,yc,yrange)

    if ( ierror == 0) then

      par(1:npar) = paropt(1:npar)
      grb(1:ncofrb) = grbopt(1:ncofrb)

      write ( *, * ) ' '
      write ( *, * ) 'Optimizing parameters:'
      call prpar(iopt,npar,nparb,nparf,par)
      write ( *, * ) ' '
      write ( *, * ) 'Optimal cost = ',cost
      write ( *, * ) ' '
      write ( *, * ) 'Number of standard full solutions:',numopt
      write ( *, * ) 'Number of auxilliary solutions:   ',numdif
    else
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  The optimization was unsuccessful.'
    end if
!
!  PAR(*) = *
!
  else if ( s_eqi ( command(1:4),'par(')) then

    call chrcti(command(5:),ival,ierror,lchar)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  ChrCTI returned nonzero error flag!'
      cycle
    end if

    if ( ival < 1 .or. maxpar < ival ) then
      write ( *, * ) ' '
      write ( *, * ) 'INPUT - Warning!'
      write ( *, * ) '  Index IVAL of PAR is out of bounds!'
      write ( *, * ) '  IVAL = ',ival
      cycle
    end if

    call chrctd(command(5+lchar+2:),value,ierror,lchar)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  ChrCTD returned nonzero error flag!'
      cycle
    end if

    par(ival) = value

    if ( 0 < iprint ) then
      write ( *, * ) 'PAR(',ival,') set to ',par(ival)
    end if
!
!  PARTAR(*) = *
!
  else if ( s_eqi ( command(1:7),'partar(')) then

    call chrcti(command(8:),ival,ierror,lchar)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  ChrCTI returned nonzero error flag!'
      cycle
    end if

    if ( ival < 1 .or. maxpar < ival ) then
      write ( *, * ) ' '
      write ( *, * ) 'INPUT - Warning!'
      write ( *, * ) '  Index IVAL of PARTAR is out of bounds!'
      write ( *, * ) '  IVAL = ', ival
      cycle
    end if

    call chrctd(command(8+lchar+2:),value,ierror,lchar)

    if ( ierror /= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  ChrCTD returned nonzero error flag!'
      cycle
    end if

    partar(ival) = value
    par(ival) = value

    if ( 0 < iprint ) then
      write ( *, * ) 'PARTAR(',ival,') set to ',partar(ival)
    end if
!
!  PICFL
!
  else if ( s_eqi ( command,'picfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - PicFL:'
      write ( *, * ) '  Apply Picard to full solution estimate GFL.'
    end if

    call picfl(afl,area,eqn,gfl,ierror,ifs,indx,ipivfl,iwrite,ldafl,maxsim, &
      nelem,neqnfl,nlband,node,np,npar,nparf,numsim,par,phifl,region,resfl, &
      rmax,splflo,tauflo,tolsim,xc,xrange,yc,yrange)

    if ( iwrite <= 1) then
      write ( *, * ) '  Picard step ',numsim,' residual norm  = ',rmax
    end if
!
!  PICRB
!
  else if ( s_eqi ( command,'picrb')) then

    if ( 0 < iprint ) then  
      write ( *, * ) 'ARBY4 - PicRB:'
      write ( *, * ) '  Apply Picard to reduced solution estimate GRB.'
    end if

    if ( ncofrb <= 0) then
      cycle
    end if

    call picrb(arb,area,grb,ierror,indx,ipivrb,iwrite,maxcofrb,maxelm,maxnfl, &
      maxsim,nbcrb,ncofrb,nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb, &
      resrb,rmax,tauflo,tolsim,xc,xrange,yc,yrange)

    if ( iwrite <= 1) then
      write ( *, * ) '  Final Picard residual was MxNorm(FXRB) = ',rmax
    end if
!
!  PRFPFL
!
  else if ( s_eqi ( command,'prfpfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr FP FL:'
      write ( *, * ) '  Print full jacobian.'
      write ( *, * ) '  Rows ILO = ',ilo,' to IHI=',ihi
      write ( *, * ) '  Cols JLO = ',jlo,' to JHI=',jhi
      write ( *, * ) ' '
      write ( *, * ) '  Parameters for matrix, PARAFL:'

      call prpar(iopt,npar,nparb,nparf,parafl)
    end if

    call prbmat(afl,ihi,ilo,jhi,jlo,ldafl,neqnfl,nlband)
!
!  PRFPRB
!
  else if ( s_eqi ( command,'prfprb')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr FP RB:'
      write ( *, * ) '  Print reduced jacobian.'
      write ( *, * ) '  Rows ILO = ',ilo,' to IHI=',ihi
      write ( *, * ) '  Cols JLO = ',jlo,' to JHI=',jhi
      write ( *, * ) ' '
      write ( *, * ) '  Parameters for matrix, PARARB:'
      call prpar(iopt,npar,nparb,nparf,pararb)
    end if

    mlo = 1
    mhi = maxcofrb
    nlo = 1
    nhi = maxcofrb

    call prdmat(arb,ihi,ilo,jhi,jlo,mhi,mlo,nhi,nlo)
!
!  PRDAT
!
  else if ( s_eqi ( command,'prdat')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr Dat'
      write ( *, * ) '  Print current problem data.'
    end if

    call prdat(disfil,drey,epsdif,gridx,gridy,hx,hy,ibs,ibump,ifs,ijac,iopt, &
      maxnew,maxopt,maxsim,nbcrb,ncofrb,nelem,nferb,neqnfl,np,npar,nparb, &
      nparf,ntay,nx,ny,region,reytay,tecfil,tolnew,tolopt,tolsim,wateb, &
      watep,wateu,watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)
!
!  PRELEM
!
  else if ( s_eqi ( command,'prelem')) then

    if ( 0 < iprint ) then  
      write ( *, * ) 'ARBY4 - Pr Elem'
      write ( *, * ) '  Print element data.'
    end if

    call prelem(ihi,ilo,nelem,node,np,xc,yc)
!
!  PRFXFL
!
  else if ( s_eqi ( command,'prfxfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr FX FL'
      write ( *, * ) '  Print full residual FXFL.'
    end if

    call prvecfl(eqn,ihi,ilo,indx,neqnfl,np,resfl)
!
!  PRFXFLNRM
!
  else if ( s_eqi ( command,'prfxflnrm')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr FX FL Nrm'
      write ( *, * ) '  Print norm of full residual FXFL.'
    end if

    call prfxfln(neqnfl,resfl)
!
!  PRFXRB
!
  else if ( s_eqi ( command,'prfxrb')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr FX RB'
      write ( *, * ) '  Print reduced residual FXRB.'
    end if

    nlo = 1
    nhi = ncofrb
    call prvecrb(ihi,ilo,nhi,nlo,resrb)
!
!  PRGFL
!
  else if ( s_eqi ( command,'prgfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr G FL:'
      write ( *, * ) '  Print full solution GFL.'
      write ( *, * ) ' '
      write ( *, * ) '  Flow parameters, PAR:'
      call prpar(iopt,npar,nparb,nparf,par)
    end if

    call prvecfl(eqn,ihi,ilo,indx,neqnfl,np,gfl)
!
!  PRGFLNRM
!
  else if ( s_eqi ( command,'prgflnrm')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr GFL Nrm:'
      write ( *, * ) '  Print norms of full solution GFL.'
    end if

    call fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar,nparf,par, &
      phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)

    call nrmflo(gfl,indx,neqnfl,np,resfl)
!
!  PRGRB
!
  else if ( s_eqi ( command,'prgrb')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr G RB:'
      write ( *, * ) '  Print reduced solution GRB.'
    end if

    call prgrb(grb,ncofrb)
!
!  PRGSEN
!
  else if ( s_eqi ( command,'prgsen')) then

    call prgrb(gsen,nbcrb+nsenfl)
!
!  PRINDX
!
  else if ( s_eqi ( command,'prindx')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr INDX'
      write ( *, * ) '  Print node/equation table,'
      write ( *, * ) '  for nodes ILO = ',ilo,' to IHI=',ihi
    end if

    call prindx(ihi,ilo,indx,np,xc,yc)
!
!  PRPAR
!
  else if ( s_eqi ( command,'prpar')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr PAR: Print current parameters PAR.'
    end if

    call prpar(iopt,npar,nparb,nparf,par)
!
!  PRPARTAR
!
  else if ( s_eqi ( command,'prpartar')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr PAR TAR:'
      write ( *, * ) '  Print target parameters PARTAR.'
    end if

    call prpar(iopt,npar,nparb,nparf,partar)
!
!  PRRBASE
!
  else if ( s_eqi ( command,'prrbase')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr RBase'
      write ( *, * ) '  Print the "R" factor of the reduced basis.'
    end if

    ilo = 1
    ihi = ncofrb
    jlo = 1
    jhi = ncofrb
    mlo = 1
    mhi = maxcofrb
    nlo = 1
    nhi = maxcofrb

    call prdmat(rbase,ihi,ilo,jhi,jlo,mhi,mlo,nhi,nlo)
!
!  PRRB
!
  else if ( s_eqi ( command,'prrb')) then

    if ( ncofrb <= 0) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  PRRB command cancelled, NCOFRB = ',ncofrb
      write ( *, * ) '  Use the GETRB command first!'
      cycle
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr RB'
      write ( *, * ) '  Print the reduced basis,'
      write ( *, * ) '  nodes ILO = ',ilo,' to IHI=',ihi
      write ( *, * ) '  columns JLO = ',jlo,' to JHI=',jhi
      write ( *, * ) ' '
      write ( *, * ) '  Parameters at reduced basis, PARRB:'
      call prpar(iopt,npar,nparb,nparf,parrb)
    end if

    call prmatfl(rb,eqn,ihi,ilo,indx,jhi,jlo,maxnfl,ncofrb,neqnfl,np)
!
!  PRSENFL: Print full sensitivities.
!
  else if ( s_eqi ( command,'prsenfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr Sen FL'
      write ( *, * ) '  Print full sensitivities.'
      write ( *, * ) ' '
      write ( *, * ) '  Parameters at sensitivity, PARSEN:'

      call prpar(iopt,npar,nparb,nparf,parsen)
    end if

    call prmatfl(senfl,eqn,ihi,ilo,indx,jhi,jlo,maxnfl,nsenfl,neqnfl,np)
!
!  PRSENNRM: Print full sensitivity norms.
!
  else if ( s_eqi ( command,'prsennrm')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr Sen Nrm'
      write ( *, * ) '  Print sensitivity norms.'
    end if

    call prsenn(maxcofrb,maxnfl,ncofrb,neqnfl,senfl)
!
!  PRSENRB: Print reduced sensitivities.
!
  else if ( s_eqi ( command,'prsenrb')) then

    if ( 0 < iprint ) then  
      write ( *, * ) 'ARBY4 - Pr Sen RB'
      write ( *, * ) '  Print matrix of reduced sensitivities,'
      write ( *, * ) '  rows ILO = ',ilo,' to IHI=',ihi
      write ( *, * ) '  columns JLO = ',jlo,' to JHI=',jhi
      write ( *, * ) ' '
      write ( *, * ) '  Parameters at sensitivity, PARSEN:'

      call prpar(iopt,npar,nparb,nparf,parsen)
    end if

    mlo = 1
    mhi = maxcofrb
    nlo = 1
    nhi = ncofrb

    call prdmat(senrb,ihi,ilo,jhi,jlo,mhi,mlo,nhi,nlo)
!
!  PRUVPGFL
!
  else if ( s_eqi ( command,'pruvpgfl')) then

    call pruvpfl(gfl,indx,neqnfl,np,xc,xmax,xmin,yc,ymax,ymin)
!
!  PRUVPRB
!
  else if ( s_eqi ( command,'pruvprb')) then

    do j = 1, ncofrb

      write ( *, * ) ' '
      write ( *, * ) 'Reduced basis vector ',j
      write ( *, * ) ' '

      call pruvpfl(rb(1,j),indx,neqnfl,np,xc,xmax,xmin,yc,ymax,ymin)

    end do
!
!  PRUVPSENFL
!
  else if ( s_eqi ( command,'pruvpsenfl')) then

    do j = 1, nsenfl

      write ( *, * ) ' '
      write ( *, * ) 'Sensitivity vector ',j
      write ( *, * ) ' '

      call pruvpfl(senfl(1,j),indx,neqnfl,np,xc,xmax,xmin,yc,ymax,ymin)

    end do
!
!  PRUVPGRB
!
  else if ( s_eqi ( command,'pruvpgrb')) then

    call pruvprb(grb,indx,maxnfl,ncofrb,nelem,node,nodelm,np, &
      rb,xc,xmax,xmin,yc,ymax,ymin)
!
!  PRXY
!
  else if ( s_eqi ( command,'prxy')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Pr XY'
      write ( *, * ) '  Print out X and Y nodal coordinates.'
    end if

    call prxy(ihi,ilo,np,ny,xc,yc)
!
!  QUIT
!
  else if ( s_eqi ( command, 'quit' ) ) then

    exit
!
!  REDUCE GFL
!
  else if ( s_eqi ( command,'reduce gfl')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Reduce GFL:'
      write ( *, * ) '  Given a reduced basis RB computed at the'
      write ( *, * ) '  full solution GFLRB, and an arbitrary full'
      write ( *, * ) '  solution GFL, compute the reduced basis '
      write ( *, * ) '  coefficients of GFL:'
      write ( *, * ) '    GRB = RB^T * GFL.'
    end if

    call gfl2rb(gfl,gflrb,grb,maxnfl,ncofrb,neqnfl,rb)
!
!  REGION =
!
  else if ( s_eqi ( command(1:6),'region')) then

    if ( s_eqi ( command(1:7),'region=')) then
      region = command(8:)
    else
      write ( *, * ) 'Enter the region, CAVITY, CHANNEL or STEP:'
      read(*,'(a)')region
      write(17,'(a)')region
    end if

    if ( s_eqi ( region,'cavity')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - Cavity:'
        write ( *, * ) '  Set user input values to cavity defaults.'
      end if

      call cavity(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb,nparf, &
        npe,nx,ny,par,region,reynld,tolnew,tolopt,tolsim,wateb,watep,wateu, &
        watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)

    else if ( s_eqi ( region,'cavity2')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - Cavity2:'
        write ( *, * ) '  Set H C Lee cavity defaults.'
      end if

      call cavity2(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb,nparf, &
        npe,nx,ny,par,region,reynld,tolnew,tolopt,tolsim,wateb,watep,wateu, &
        watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)

    else if ( s_eqi ( region,'channel')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - Channel:'
        write ( *, * ) '  Set user input values to channel defaults.'
      end if

      call channl(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb,nparf, &
        npe,nx,ny,par,region,reynld,tolnew,tolopt,tolsim,wateb,watep,wateu, &
        watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)

    else if ( s_eqi ( region,'step')) then

      if ( 0 < iprint ) then
        write ( *, * ) 'ARBY4 - Step`:'
        write ( *, * ) '  Set user input values to step defaults.'
      end if

      call step(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb,nparf,npe, &
        nx,ny,par,region,reynld,tolnew,tolopt,tolsim,wateb,watep,wateu,watev, &
        xbl,xbr,xprof,xrange,ybl,ybr,yrange)

    end if
!
!  REYNLD =
!
  else if ( s_eqi ( command(1:6),'reynld')) then

    if ( s_eqi ( command(1:7),'reynld=')) then
      read(command(8:),*)reynld
    else
      write ( *, * ) 'Enter value for REYNLD:'
      read(*,*)reynld
      write(17,*)reynld
    end if

    par(nparf+nparb+1) = reynld

    if ( 0 < iprint ) then
      write ( *, * ) 'REYNLD parameter set to ',reynld
    end if
!
!  REYTAY =
!
  else if ( s_eqi ( command(1:6),'reytay')) then

    if ( s_eqi ( command(1:7),'reytay=')) then
      if ( s_eqi ( command,'reytay=reynld')) then
        reytay = reynld
      else
        read(command(8:),*)reytay
      end if
    else
      write ( *, * ) 'Enter value for REYTAY:'
      read(*,*)reytay
      write(17,*)reytay
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'REYTAY parameter set to ',reytay
    end if
!
!  SETGEO
!
  else if ( s_eqi ( command,'setgeo')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - SetGeo: Set problem geometry.'
    end if

    call setgeo(area,etaq,gridx,gridy,ibs,isotri,nelem,node,nodelm,np,npar, &
      nparb,nparf,nx,ny,par,phifl,region,splbmp,taubmp,wquad,xbl,xbr,xc, &
      xquad,xrange,xsiq,ybl,ybr,yc,yquad,yrange)

!
!  SETLOG
!
  else if ( s_eqi ( command,'setlog')) then

    if ( 0 < iprint ) then  
      write ( *, * ) 'ARBY4 - SetLog: Set problem logical data.'
    end if

    call setlog(eqn,hx,hy,ibump,indx,isotri,ldafl,maxelm,maxnfl,maxnp,nelem, &
      neqnfl,nlband,node,np,nprof,nx,ny,region,xbl,xbr,xprof,xrange,ybr,yrange)
!
!  STOP
!
  else if ( s_eqi ( command,'stop') ) then

    exit
!
!  TARGET
!
  else if ( s_eqi ( command,'target')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Target:'
      write ( *, * ) '  Save current GFL as GTAR.'
    end if

    call target(cost0,gfl,gfltar,indx,maxnfl,maxny,maxparb,neqnfl,np,npar, &
      nparb,nprof,ny,par,partar,splbmp,taubmp,wateb,watep,wateu,watev,xbl, &
      xbr,ybl,ybr,yc)
!
!  TEST2
!
  else if ( s_eqi ( command,'test2')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Test2:'
      write ( *, * ) '  Compare full and reduced state variables'
      write ( *, * ) '  in elements ILO through IHI.'
    end if

    call test2(gfl,grb,ihi,ilo,indx,maxcofrb,maxelm,ncofrb, &
      nelem,neqnfl,node,np,phifl,phirb)
!
!  TEST3
!
  else if ( s_eqi ( command,'test3')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Test3:'
      write ( *, * ) '  Compare RB*Rfact and SenFL.'
    end if

    call test3(maxcofrb,maxnfl,ncofrb,neqnfl,rb,senfl,senrb)
!
!  TEST4
!
  else if ( s_eqi ( command,'test4')) then

    if ( ipar <= 0 .or. npar < ipar ) then
      write ( *, * ) ' '
      write ( *, * ) 'ARBY4 - Warning!'
      write ( *, * ) '  Cancelling the TEST4 command.'
      write ( *, * ) '  IPAR = ',ipar
      write ( *, * ) '  but IPAR must be at least 0'
      write ( *, * ) '  and no more than NPAR = ',npar
      cycle
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Test4'
      write ( *, * ) '  Compare full sensitivities computed directly'
      write ( *, * ) '  and via finite differences.'
    end if

    call test4(afl,area,difcof,dpar,drey,eqn,gfl,gflafl, &
      ifs,ijac,indx,ipar,ipivfl,iwrite,ldafl,maxcofrb,maxelm, &
      maxnew,maxnfl,ncofrb,nelem,neqnfl,nlband,node,np,npar, &
      nparf,nsenfl,par,parafl,phifl,region,resfl, &
      senfl,splflo,tauflo,tolnew,xrange,yc,yrange)
!
!  TEST5
!
  else if ( s_eqi ( command,'test5')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Test5:'
      write ( *, * ) '  Compare RB*Rfact and basis vectors.'
    end if

    call test5(maxcofrb,maxnfl,ncofrb,neqnfl,rb,rbase)
!
!  TIME
!
  else if ( s_eqi ( command,'time')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - Time:'
      write ( *, * ) '  Report current time,'
      write ( *, * ) '  time elapsed since last TIME call,'
      write ( *, * ) '  time elapsed since the program began.'
    end if

    write ( *, * ) '  The (real) start time was  '// trim ( tstart )
    call date_and_time ( date, time )
    write ( *, * ) '  The current (real) time is ' // trim ( time )
!
!  TOLNEW =
!
  else if ( s_eqi ( command(1:6),'tolnew')) then

    if ( s_eqi ( command(1:7),'tolnew=')) then
      read(command(8:),*)tolnew
    else
      write ( *, * ) 'Enter value for TOLNEW:'
      read(*,*)tolnew
      write(17,*)tolnew
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'TOLNEW set to ',tolnew
    end if
!
!  TOLOPT =
!
  else if ( s_eqi ( command(1:6),'tolopt')) then

    if ( s_eqi ( command(1:7),'tolopt=')) then
      read(command(8:),*)tolopt
    else
      write ( *, * ) 'Enter value for TOLOPT:'
      read(*,*)tolopt
      write(17,*)tolopt
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'TOLOPT set to ',tolopt
    end if
!
!  TOLSIM =
!
  else if ( s_eqi ( command(1:6),'tolsim')) then

    if ( s_eqi ( command(1:7),'tolsim=')) then
      read(command(8:),*)tolsim
    else
      write ( *, * ) 'Enter value for TOLSIM:'
      read(*,*)tolsim
      write(17,*)tolsim
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'TOLSIM set to ',tolsim
    end if
!
!  TecPlot
!
  else if ( s_eqi ( command,'tecplot')) then

    if ( 0 < iprint ) then
      write ( *, * ) 'ARBY4 - TecPlot:'
      write ( *, * ) '  Write data to TECPLOT plot file.'
    end if

    call intprs(gfl,indx,nelem,neqnfl,node,np,p)

    u(1:np) = gfl(indx(1,1:np))
    v(1:np) = gfl(indx(2,1:np))

    call wrtec(nelem,node,np,p,tecfil,u,v,xc,yc)
!
!  WATEB =
!
  else if ( s_eqi ( command(1:5),'wateb')) then

    if ( s_eqi ( command(1:6),'wateb=')) then
      read(command(7:),*)wateb
    else
      write ( *, * ) 'Enter value for WATEB:'
      read(*,*)wateb
      write(17,*)wateb
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'WATEB set to ',wateb
    end if
!
!  WATEP =
!
  else if ( s_eqi ( command(1:5),'watep')) then

    if ( s_eqi ( command(1:6),'watep=')) then
      read(command(7:),*)watep
    else
      write ( *, * ) 'Enter value for WATEP:'
      read(*,*)watep
      write(17,*)watep
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'WATEP set to ',watep
    end if
!
!  WATEU =
!
  else if ( s_eqi ( command(1:5),'wateu')) then

    if ( s_eqi ( command(1:6),'wateu=')) then
      read(command(7:),*)wateu
    else
      write ( *, * ) 'Enter value for WATEU:'
      read(*,*)wateu
      write(17,*)wateu
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'WATEU set to ',wateu
    end if
!
!  WATEV =
!
  else if ( s_eqi ( command(1:5),'watev')) then

    if ( s_eqi ( command(1:6),'watev=')) then
      read(command(7:),*)watev
    else
      write ( *, * ) 'Enter value for WATEV:'
      read(*,*)watev
      write(17,*)watev
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'WATEV set to ',watev
    end if
!
!  XBL =
!
  else if ( s_eqi ( command(1:3),'xbl')) then

    if ( s_eqi ( command(1:4),'xbl=')) then
      read(command(5:),*)xbl
    else
      write ( *, * ) 'Enter value for XBL:'
      read(*,*)xbl
      write(17,*)xbl
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'XBL set to ',xbl
    end if
!
!  XBR =
!
  else if ( s_eqi ( command(1:3),'xbr')) then

    if ( s_eqi ( command(1:4),'xbr=')) then
      read(command(5:),*)xbr
    else
      write ( *, * ) 'Enter value for XBR:'
      read(*,*)xbr
      write(17,*)xbr
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'XBR set to ',xbr
    end if
!
!  XMAX =
!
  else if ( s_eqi ( command(1:4),'xmax')) then

    if ( s_eqi ( command(1:5),'xmax=')) then
      read(command(6:),*)xmax
    else
      write ( *, * ) 'Enter value for XMAX:'
      read(*,*)xmax
      write(17,*)xmax
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'XMAX set to ',xmax
    end if
!
!  XMIN =
!
  else if ( s_eqi ( command(1:4),'xmin')) then

    if ( s_eqi ( command(1:5),'xmin=')) then
      read(command(6:),*)xmin
    else
      write ( *, * ) 'Enter value for XMIN:'
      read(*,*)xmin
      write(17,*)xmin
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'XMIN set to ',xmin
    end if
!
!  XPROF =
!
  else if ( s_eqi ( command(1:5),'xprof')) then

    if ( s_eqi ( command(1:6),'xprof=')) then
      read(command(7:),*)xprof
    else
      write ( *, * ) 'Enter value for XPROF:'
      read(*,*)xprof
      write(17,*)xprof
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'XPROF set to ',xprof
    end if
!
!  XRANGE =
!
  else if ( s_eqi ( command(1:6),'xrange')) then

    if ( s_eqi ( command(1:7),'xrange=')) then
      read(command(8:),*)xrange
    else
      write ( *, * ) 'Enter value for XRANGE:'
      read(*,*)xrange
      write(17,*)xrange
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'XRANGE set to ',xrange
    end if
!
!  YBL =
!
  else if ( s_eqi ( command(1:3),'ybl')) then

    if ( s_eqi ( command(1:4),'ybl=')) then
      read(command(5:),*)ybl
    else
      write ( *, * ) 'Enter value for YBL:'
      read(*,*)ybl
      write(17,*)ybl
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'YBL set to ',ybl
    end if
!
!  YBR =
!
  else if ( s_eqi ( command(1:3),'ybr')) then

    if ( s_eqi ( command(1:4),'ybr=')) then
      read(command(5:),*)ybr
    else
      write ( *, * ) 'Enter value for YBR:'
      read(*,*)ybr
      write(17,*)ybr
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'YBR set to ',ybr
    end if
!
!  YMAX =
!
  else if ( s_eqi ( command(1:4),'ymax')) then

    if ( s_eqi ( command(1:5),'ymax=')) then
      read(command(6:),*)ymax
    else
      write ( *, * ) 'Enter value for YMAX:'
      read(*,*)ymax
      write(17,*)ymax
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'YMAX set to ',ymax
    end if
!
!  YMIN =
!
  else if ( s_eqi ( command(1:4),'ymin')) then

    if ( s_eqi ( command(1:5),'ymin=')) then
      read(command(6:),*)ymin
    else
      write ( *, * ) 'Enter value for YMIN:'
      read(*,*)ymin
      write(17,*)ymin
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'YMIN set to ',ymin
    end if
!
!  YRANGE =
!
  else if ( s_eqi ( command(1:6),'yrange')) then

    if ( s_eqi ( command(1:7),'yrange=')) then
      read(command(8:),*)yrange
    else
      write ( *, * ) 'Enter value for YRANGE:'
      read(*,*)yrange
      write(17,*)yrange
    end if

    if ( 0 < iprint ) then
      write ( *, * ) 'YRANGE set to ',yrange
    end if
!
!  Unrecognized command
!
  else

    write ( *, * ) ' '
    write ( *, * ) 'ARBY4 - Warning!'
    write ( *, * ) '  Unrecognized command: ' // trim ( command )

  end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'ARBY4 - STOP command:'
  write ( *, * ) '  Halt the program!'
  write ( *, * ) ' '

  close(unit = 17)

  write ( *, * ) '  Closing the user input file ARBY.IN.'
  write ( *, * ) ' '
  write ( *, * ) '  The (real) start time was    '// trim ( tstart )

  call date_and_time ( date, time )
  tstop = time

  write ( *, * ) '  The (real) stopping time was '// trim ( tstop )

  call delhms ( tstart, tstop, itemp )

  write ( *, * ) '  The (real) elapsed time in seconds is ',itemp
  write ( *, * ) '  The real elapsed time in minutes is ', &
    real(itemp) / 60.0D+00
  call etime ( tarray, estop )
  write ( *, * ) ' '
  write ( *, * ) '  CPU in seconds = ',estop-estart
  write ( *, * ) '  CPU in minutes = ',(estop-estart)/60.0D+00
!
!  Terminate.
!
  write ( *, * ) ' '
  write ( *, * ) 'ARBY4:'
  write ( *, * ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine diffprb ( arb, area, epsdif, grb, indx, maxcofrb, maxelm, &
  maxnfl, nbcrb, ncofrb, nelem, nferb, node, np, npar, nparf, nx, &
  ny, par, phirb, rb, resrb, reynld, tauflo, xc, xrange, yc, yrange )

!*****************************************************************************80
!
!! DIFFPRB estimates the jacobian of the reduced function.
!
!  Discussion:
!
!    Finite differences are used.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) ARB(MAXCOFRB,MAXCOFRB).
!    ARB contains the Jacobian or Picard matrix for the reduced
!    Navier Stokes system, stored as an NCOFRB by NCOFRB array.
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!
!    or, if the element is isoperimetric,
!
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    real ( kind = 8 ) EPSDIF.
!    EPSDIF is a small quantity, which is used to compute the
!    perturbations for the finite difference approximations.
!
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    If K = INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry GFL(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either because
!    the variable is specified in some other way, or because
!    (in the case of pressure), there is no coefficient associated
!    with that node.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXELM.
!    MAXELM is the maximum number of elements.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NFERB.
!    NFERB is the number of reduced basis coefficients that will
!    be determined via the finite element method.
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    The local ordering of the nodes is suggested by this diagram:
!
!      Global nodes   Elements      NODE
!                                                     1  2  3  4  5  6
!      74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                     |  /   /|
!      73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                     |/   /  |
!      72  82  92     2   1-6-3
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!    The number of elements along a line in the X direction is
!    NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!    integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!    The number of elements along a line in the Y direction is
!    NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    real ( kind = 8 ) PHIRB(3,MAXCOFRB,15,MAXELM).
!    PHIRB contains the values of a finite element basis function
!    or its X or Y derivative, in a given element, at a given
!    quadrature point, for a particular reduced basis function.
!    For PHIRB(I,J,K,L), index J refers to the reduced basis
!    basis functions, for J = 0 to NCOFRB.
!    The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!      For the quadrature point I, and reduced basis function J,
!      in element L, PHIRB(I,J,K,L) represents the value of:
!
!        K = 1, WUrb, the finite element U velocity basis function;
!        K = 2, dWUrbdX, the X derivative of WUrb;
!        K = 3, dWUrbdY, the Y derivative of WUrb;
!        K = 4, WVrb, the finite element V velocity basis function;
!        K = 5, dWVrbdX, the X derivative of WVrb;
!        K = 6, dWVrbdY, the Y derivative of WVrb;
!        K = 7, Q, the finite element pressure basis function.
!        K = 8, dQrbdX, the X derivative of Qrb;
!        K = 9, dQrbdY, the Y derivative of Qrb.
!        K = 10, WU0rb, same as WUrb, with zero BC.
!        K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!        K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!        K = 13, WV0rb, same as WVrb, with zero BC.
!        K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!        K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!    RB is based on a particular solution of the full system,
!    which is saved as GFLRB.
!    We compute the first NCOFRB derivatives of GFLRB with
!    respect to a parameter.  The first derivative
!    is stored in column 1 of RB, and so on. 
!
!    real ( kind = 8 ) RESRB(NCOFRB).
!    RESRB contains the residual in the reduced basis equations,
!    for the parameter values PAR and reduced basis coefficients GRB.
!
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the current value of the Reynolds number.
!    Normally, REYNLD is stored as PARA(NPARF+NPARB+1).
!
!    real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow. 
!
!    A recent code change was made.  For a channel flow, where
!    NPARF = 1 meant a fit through 3 points, now NPARF=3 means
!    a fit through 3 points.  The endpoints must be explicitly
!    counted.
!
!    real ( kind = 8 ) XC(NP).
!    XC contains the X coordinates of the nodes.
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none

  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf

  real ( kind = 8 ) arb(maxcofrb,ncofrb)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) delta
  real ( kind = 8 ) epsdif
  real ( kind = 8 ) grb(ncofrb)
  real ( kind = 8 ) grbtmp(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,maxelm)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Get the function value at the base value.
!
  call fxrb(area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
    nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb, &
    resrb,reynld,tauflo,xc,xrange,yc,yrange)
!
!  Start each column of the jacobian equal to F(GRB).
!
  do j = 1, ncofrb
    arb(1:ncofrb,j) = resrb(1:ncofrb)
  end do

  do j = 1, ncofrb
!
!  Perturb G(J).
!
    grbtmp(1:ncofrb) = grb(1:ncofrb)
    delta = epsdif * ( 1.0D+00 + abs ( grb(j) ) )
    grbtmp(j) = grb(j) + delta
!
!  Evaluate F(I) at the perturbed value G(J).
!
    call fxrb(area,grbtmp,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
      nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb, &
      resrb,reynld,tauflo,xc,xrange,yc,yrange)
!
!  Estimate the dependence ARB(I,J) = dF(I)/dG(J)
!
    arb(1:ncofrb,j) = ( resrb(1:ncofrb) - arb(1:ncofrb,j) ) / delta

  end do

  return
end
subroutine difsenfl ( afl, area, difcof, dpar, eqn, gfl, gflafl, ifs, &
  ijac, indx, ipar, ipivfl, iwrite, ldafl, maxcofrb, maxelm, maxnew, &
  maxnfl, ncofrb, nelem, neqnfl, nlband, node, np, npar, nparf, par, &
  parafl, phifl, region, resfl, senfl, splflo, tauflo, tolnew, &
  xrange, yc, yrange )

!*****************************************************************************80
! 
!! DIFSENFL estimates full solution derivatives with respect to parameters.
!
!  Discussion:
!
!    The routine computes a central difference estimate for the first 
!    NCOFRB derivatives of the full solution GFL with respect to the IPAR-th
!    parameter.
!
!    DIFSENFL is rather inefficient.  ALTHOUGH SOME SOLUTIONS
!    ARE USED SEVERAL TIMES, DIFSENFL RECOMPUTES THEM EACH TIME.
!    A CORRECTION OF THIS PROBLEM WOULD BE TO COMPUTE THE ENTIRE
!    TRIANGLE OF COEFFICIENTS FIRST, AND THEN COMPUTE JUST THE
!    SOLUTIONS NEEDED ONCE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
! 
!    real ( kind = 8 ) AFL(LDAFL,MAXNFL).
!    If Newton iteration is being carried out, AFL contains the
!    Jacobian matrix for the full system.
!    If Picard iteration is being carried out, AFL contains the
!    Picard matrix for the full system.
!    AFL is stored in LINPACK general band storage mode, with
!    logical dimensions (3*NLBAND+1, NEQNFL).
!    Where is the (I,J) entry of AFL actually stored?
!    AFL has actual storage for such an entry only if
!      -NLBAND <= I-J <= NLBAND.
!    In such a case, the (I,J) entry is actually stored in
!      AFL(I-J+2*NLBAND+1,J)
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    real ( kind = 8 ) DIFCOF(NDIF).
!    DIFCOF contains the coefficients needed to approximate
!    the 0-th through (NDIF-1)-th derivatives of a function F.
!
!    real ( kind = 8 ) DREY.
!    DREY is the suggested increment in the REYNLD value,
!    to be used during the finite difference estimations.
!
!    real ( kind = 8 ) DOPT(NPAR).
!    DOPT contains a set of scale factors for the parameters, used
!    by the optimization code.  The suggestion is that DOPT(I) be
!    chosen so that DOPT(I)*PAR(I) is roughly the same order of
!    magnitude for I from 1 to NPAR.
!
!    real ( kind = 8 ) DPAR.
!    DPAR is the suggested increment in the parameter value,
!    to be used during the finite difference estimations.
!
!    character ( len = 2 ) EQN(MAXNFL).
!    EQN records the "type" of each equation that will be generated, and
!    which is associated with an unknown. 
!
!    'U'  A horizontal momentum equation.
!    'UB' The condition U = 0 applied at a node on the bump.
!    'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U = 0 applied at a node on a fixed wall.
!    'U0' A dummy value of U = 0 should be set.
!
!    'V'  A vertical momentum equation.
!    'VB' The condition V = 0 applied at a node on the bump.
!    'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V = 0 applied at a node on a fixed wall.
!    'V0' A dummy value of V = 0 should be set.
!
!    'P'  A continuity equation.
!    'PB' The condition P = 0 applied at (XMAX,YMAX).
!    'P0' A dummy value of P = 0 should be set.
!
!    real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
!
!    real ( kind = 8 ) GFLAFL(NEQNFL).
!    GFLAFL stores the value of GFL at which the Jacobian
!    was generated.
!
!    integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) IJAC.
!    IJAC determines the frequency for evaluating and factoring
!    the Jacobian matrix during any particular Newton process.
!    1, evaluate the Jacobian on every step of the Newton
!       iteration.
!    n, evaluate the Jacobian only at steps 0, n, 2*n, and so on.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    If K = INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry GFL(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either because
!    the variable is specified in some other way, or because
!    (in the case of pressure), there is no coefficient associated
!    with that node.
!
!    integer ( kind = 4 ) IPAR.
!    IPAR is the index of the parameter to be varied.
!
!    integer ( kind = 4 ) IPIVFL(NEQNFL).
!    IPIVFL is a pivot vector for the solution of the full
!    linear system.
!
!    integer ( kind = 4 ) IWRITE.
!    IWRITE controls the amount of output printed.
!    0, print out the least amount.
!    1, print out some.
!    2, print out a lot.
!
!    integer ( kind = 4 ) LDAFL.
!    LDAFL is the first dimension of the matrix AFL as declared in
!    the main program.  LDAFL must be at least 3*NLBAND+1.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXELM.
!    MAXELM is the maximum number of elements.
!
!    integer ( kind = 4 ) MAXNEW.
!    MAXNEW is the maximum number of steps to take in one Newton
!    iteration.  A typical value is 20.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NLBAND.
!    NLBAND is the lower bandwidth of the matrix AFL.
!    The zero structure of AFL is assumed to be symmetric, and so
!    NLBAND is also the upper bandwidth of AFL.
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    The local ordering of the nodes is suggested by this diagram:
!
!      Global nodes   Elements      NODE
!                                                     1  2  3  4  5  6
!      74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                     |  /   /|
!      73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                     |/   /  |
!      72  82  92     2   1-6-3
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    real ( kind = 8 ) PARAFL(NPAR).
!    PARAFL contains the parameters where the Picard matrix or
!    Jacobian of the full system was generated.
!
!    real ( kind = 8 ) PHIFL(3,6,10,NELEM).
!    PHIFL contains the value of a finite element basis function, its
!    derivative, or other information, evaluated at the quadrature
!    points (which are the element midside nodes).
!
!    The meaning of the entry PHIFL(I,J,K,L) is as follows.
!    For the quadrature point I, and basis function J, in element L,
!    PHIFL(I,J,K,L) represents the value of:
!
!      K =  1, W, the finite element basis function for velocities;
!      K =  2, dWdX, the X derivative of W;
!      K =  3, dWdY, the Y derivative of W;
!      K =  4, Q, the finite element basis function for pressures;
!      K =  5, dQdX, the X derivative of Q;
!      K =  6, dQdY, the Y derivative of Q;
!      K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!      K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!      K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!      K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!    In particular, PHIFL(I,J,K,L) is the value of the quadratic
!    basis function W associated with local node J in element L,
!    evaluated at quadrature point I.
!
!    Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!    since there are only three linear basis functions.
!
!    character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottome, with tangential velocity specifications
!    there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    real ( kind = 8 ) RESFL(NEQNFL).
!    RESFL contains the residual in the full basis equations.
!
!    real ( kind = 8 ) SENFL(MAXNFL,MAXCOFRB).
!    Columns 1 through NSENFL of SENFL contain the sensitivities
!    of the full solution with respect to the REYNLD parameter, for
!    orders 0 through NSENFL-1.
!
!    SENFL(I,J) contains the (J-1)-th sensitivity of the I-th full unknown
!    with respect to REYNLD.
!
!    real ( kind = 8 ) SPLFLO(NPARF).
!    SPLFLO contains the spline coefficients for the inflow.
!
!    real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow. 
!
!    A recent code change was made.  For a channel flow, where
!    NPARF = 1 meant a fit through 3 points, now NPARF=3 means
!    a fit through 3 points.  The endpoints must be explicitly
!    counted.
!
!    real ( kind = 8 ) TOLNEW.
!    TOLNEW is the convergence tolerance for the Newton iteration.
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) afl(ldafl,neqnfl)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) difcof(ncofrb)
  real ( kind = 8 ) dpar
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gflafl(neqnfl)
  real ( kind = 8 ) gfltmp(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivfl(neqnfl)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) ndif
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) numnew
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) parafl(npar)
  real ( kind = 8 ) partmp(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Zero out the SENFL array.
!
  senfl(1:neqnfl,1:ncofrb) = 0.0D+00

  write ( *, * ) ' '
  write ( *, * ) '  DIFSENFL: DPAR = ',dpar
!
!  Compute difference NDIF, for NDIF = 0 to NCOFRB.
!
  do ndif = 1, ncofrb

    if ( ndif == 0 ) then

      senfl(1:neqnfl,ndif) = gfl(1:neqnfl)

    else

      if ( 2 <= iwrite ) then
        write ( *, '(a)' ) ' '
        write ( *, * ) 'DIFSENFL - Computing difference NDIF = ',ndif
      end if
!
!  Get the NDIF-1 order difference coefficients.
!
      call difset(difcof,dpar,iwrite,ndif)
!
!  Evaluate the solution at several values of the parameter.
!
      do i = 0, ndif
!
!  Copy the parameters, but reset the IPAR-th parameter value.
!
        partmp(1:npar) = par(1:npar)
        partmp(ipar) = par(ipar)+(2*i-ndif)*dpar
        write ( *, * ) 'J = ',j,' PAR(IPAR)=',partmp(ipar)
!
!  Estimate the solution GTMP at parameters PARTMP.
!
        gfltmp(1:neqnfl) = gfl(1:neqnfl)
!
!  Call NEWTFL to get the solution more closely.
!
        call newtfl(afl,area,eqn,gfltmp,gflafl,ierror,ifs,ijac,indx,ipivfl, &
          iwrite,ldafl,maxelm,maxnew,nelem,neqnfl,nlband,node,np,npar,nparf, &
          numnew,partmp,parafl,phifl,region,resfl,rmax,splflo,tauflo,tolnew, &
          xrange,yc,yrange)

        if ( ierror /= 0) then
          write ( *, * ) ' '
          write ( *, * ) 'DIFSENFL - Fatal error!'
          write ( *, * ) '  NEWTFL failed, with IERROR = ',ierror
          stop
        end if
!
!  Add the term associated with this solution to the estimate
!  of the NDIF-th derivative.
!
        senfl(1:neqnfl,ndif) = senfl(1:neqnfl,ndif) &
          + difcof(i+1) * gfltmp(1:neqnfl)

      end do
    end if
  end do

  return
end
subroutine difsenrb(arb,area,difcof,dpar,grb,grbarb,indx,ipar,ipivrb,iwrite, &
  maxcofrb,maxelm,maxnew,maxnfl,nbcrb,ncofrb,nelem,nferb,node,np,npar,nparf, &
  nx,ny,par,pararb,phirb,rb,resrb,senrb,tauflo,tolnew,xc,xrange,yc,yrange)

!*****************************************************************************80
! 
!! DIFSENRB estimates the reduced sensitivities using finite differences.
!
!  Discussion:
!
!    The denominators used in the difference calculations
!    get very small if the original increment is smaller than 1.
!    It is strongly suggested that the parameter increment not be
!    made too small!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 8 ) ARB(MAXNRB,MAXNRB).
!    ARB contains the Jacobian or Picard matrix for the reduced
!    Navier Stokes system, stored as a dense NCOFRB by NCOFRB
!    array.
!
!    Input, real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Workspace, real ( kind = 8 ) DCOF(0:NDIF).
!    DCOF contains the coefficients needed to approximate
!    the NDIF-th derivative of a function F.
!
!    Input, real ( kind = 8 ) DPAR.
!    DPAR is the suggested increment in the parameter value,
!    to be used during the finite difference estimations.
!
!    Input, real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    Output, real ( kind = 8 ) GRBARB(NCOFRB).
!    GRBARB contains the reduced basis coefficients at which
!    the matrix ARB was last evaluated.
!
!    Workspace, real ( kind = 8 ) GRBTMP(NCOFRB).
!
!    Input, integer ( kind = 4 ) IPAR.
!    The index of the parameter to be varied.
!
!    Workspace, integer IPIVRB(NCOFRB).
!    IPIVRB is a pivot vector for the solution of the reduced
!    linear system.
!
!    Input, integer ( kind = 4 ) IWRITE.
!    IWRITE controls the amount of output printed.
!    0, print out the least amount.
!    1, print out some.
!    2, print out a lot.
!
!    Input, integer ( kind = 4 ) MAXNEW.
!    MAXNEW is the maximum number of steps to take in one Newton
!    iteration.  A typical value is 20.
!
!    Input, integer ( kind = 4 ) MAXNRB.
!    MAXNRB is the maximum number of equations allowed for the
!    reduced basis system.
!
!    Input, integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    Input, integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of basis functions used for the
!    reduced basis method.  (The first basis vector is labeled
!    "0").  In this program, that amounts to the number of columns
!    in the matrix RB.  NCOFRB is also the number of reduced basis
!    state equations, and reduced basis coefficients GRB.
!
!    Input, integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!    NPAR = NPARF + NPARB + 1.
!    The parameters control the shape of the inflow,
!    the shape of the bump obstacle, and the strength of the
!    flow.
!
!    Input, real ( kind = 8 ) PAR(NPAR).
!    PAR is the current estimate for the parameters.
!    PAR(1:NPARF)             = inflow controls.
!    PAR(NPARF+1:NPARF+NPARB) = bump controls.
!    PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    Output, real ( kind = 8 ) PARARB(NPAR).
!    PARARB contains the parameters where the Picard matrix or
!    Jacobian of the reduced system was generated.
!
!    Workspace, real ( kind = 8 ) PARTMP(NPAR).
!
!    real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!    PHIRB contains the values of a finite element basis function
!    or its X or Y derivative, in a given element, at a given
!    quadrature point, for a particular reduced basis function.
!
!    For PHIRB(I,J,K,L), index J refers to the reduced basis
!    basis functions, for J = 0 to NCOFRB.
!
!    The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!      For the quadrature point I, and reduced basis function J,
!      in element L, PHIRB(I,J,K,L) represents the value of:
!
!        K = 1, WUrb, the finite element U velocity basis function;
!        K = 2, dWUrbdX, the X derivative of WUrb;
!        K = 3, dWUrbdY, the Y derivative of WUrb;
!        K = 4, WVrb, the finite element V velocity basis function;
!        K = 5, dWVrbdX, the X derivative of WVrb;
!        K = 6, dWVrbdY, the Y derivative of WVrb;
!        K = 7, Q, the finite element pressure basis function.
!        K = 8, dQrbdX, the X derivative of Qrb;
!        K = 9, dQrbdY, the Y derivative of Qrb.
!        K = 10, WU0rb, same as WUrb, with zero BC.
!        K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!        K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!        K = 13, WV0rb, same as WVrb, with zero BC.
!        K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!        K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    Workspace, real ( kind = 8 ) RESRB(NCOFRB).
!    RESRB contains the residual in the reduced basis equations,
!    for the parameter values PAR and reduced basis coefficients GRB.
!
!    Output, real ( kind = 8 ) SENRB(MAXNRB,NCOFRB).
!    SENRB contains the first several order sensitivities of the
!    reduced solution with respect to the REYNLD parameter.
!    SENRB(I,J) contains the J-th sensitivity of the I-th reduced unknown
!    with respect to REYNLD.
!
!    Input, real ( kind = 8 ) TOLNEW.
!    TOLNEW is the convergence tolerance for the Newton iteration.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) difcof(maxcofrb)
  real ( kind = 8 ) dpar
  real ( kind = 8 ) grb(maxcofrb)
  real ( kind = 8 ) grbarb(maxcofrb)
  real ( kind = 8 ) grbtmp(maxcofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icof
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivrb(maxcofrb)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jcof
  integer ( kind = 4 ) jdif
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) pararb(npar)
  real ( kind = 8 ) partmp(npar)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) resrb(maxcofrb)
  real ( kind = 8 ) senrb(maxcofrb,maxcofrb)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Zero out the SENRB array.
!
  senrb(1:maxcofrb,1:maxcofrb) = 0.0D+00
!
!  JCOF counts the number of coefficients we will compute on
!  each pass.  We're done on the last pass.
!
  do jcof = 1, ncofrb

    jdif = jcof-1

    if ( jdif == 0) then

      senrb(1:ncofrb,jcof) = grb(1:ncofrb)

    else
      write ( *, * ) ' '
      write ( *, * ) 'Computing difference order JDIF = ',jdif
!
!  Get the JCOF difference coefficients DIFCOF.
!
      call difset(difcof,dpar,iwrite,jcof)
!
!  Evaluate the solution at JCOF values of the parameter.
!
      do icof = 1, jcof
!
!  Copy the parameters, but reset the IPAR-th parameter value.
!
        partmp(1:npar) = par(1:npar)
        partmp(ipar) = par(ipar)+(2*icof-jcof-1)*dpar
        write ( *, * ) 'ICOF = ',ICOF,' PAR(IPAR)=',partmp(ipar)
!
!  Estimate the solution GRBTMP at parameters PARTMP.
!
        grbtmp(1:ncofrb) = grb(1:ncofrb)
!
!  Call NEWTRB to get the solution more closely.
!
        write ( *, * ) 'About to call NEWTRB'

        call newtrb(arb,area,grbtmp,grbarb,ierror,indx,ipivrb,iwrite, &
          maxcofrb,maxelm,maxnew,maxnfl,nbcrb,ncofrb,nelem,nferb,node,np, &
          npar,nparf,nx,ny,partmp,pararb,phirb,rb,resrb,rmax,tauflo,tolnew, &
          xc,xrange,yc,yrange)

        if ( ierror /= 0) then
          write ( *, * ) ' '
          write ( *, * ) 'DIFSENRB - Fatal error!'
          write ( *, * ) '  NEWTRB failed, with IERROR = ',ierror
          stop
        end if
!
!  Add the term associated with this solution to the estimate
!  of the JDIF-th derivative.
!
        write ( *, * ) 'ABOUT TO ADD BLEEDING TERM'
        do j = 1, ncofrb
          senrb(j,jcof) = senrb(j,jcof)+difcof(icof)*grbtmp(j)
        end do

      end do

    end if

  end do

  return
end
subroutine flowbc(ifs,npar,nparf,par,region,splflo,tauflo,ubc,vbc, &
  xrange,xval,yrange,yval)

!*****************************************************************************80
!
!! FLOWBC computes the specified boundary values at a given position.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    Input, integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow, the
!    shape of the bump, and the value of the Reynolds number.
!
!    Input, integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    Input, real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    Input, character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottom, with tangential velocity specifications
!    there.
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    Workspace, real ( kind = 8 ) SPLFLO(NPARF).
!    SPLFLO contains the spline coefficients for the inflow.
!
!    Workspace, real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow.  There are NPARF of them, because the end
!    values of the spline are constrained to have particular
!    values.
!
!    Output, real ( kind = 8 ) UBC.
!    UBC is the value of the horizontal velocity specified
!    at (XVAL,YVAL).
!
!    Output, real ( kind = 8 ) VBC.
!    VBC is the value of the vertical velocity specified at (XVAL,YVAL).
!
!    Input, real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    Input, real ( kind = 8 ) XVAL.
!    XVAL is the X coordinate of the point on the inflow boundary at
!    which the specified velocity is desired.
!
!    Input, real ( kind = 8 ) YRANGE.
!    YRANGE is the total width of the region.
!
!    Input, real ( kind = 8 ) YVAL.
!    YVAL is the Y coordinate of the point on the inflow boundary at
!    which the specified velocity is desired.
!
  implicit none
!
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifs
  logical s_eqi
  real ( kind = 8 ) par(npar)
  character ( len = 20 ) region
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xval
  real ( kind = 8 ) yrange
  real ( kind = 8 ) yval
!
!  Inflow points for the cavity have the form (X,YRANGE).
!  NPARF must be at least 1.
!
  if ( s_eqi ( region,'cavity')) then

    if ( nparf == 1) then
      tauflo(1) = xrange/2.0D+00
    else
      do i = 1, nparf
        tauflo(i) = xrange*real ( i - 1, kind = 8 )/ real ( nparf-1, kind = 8 )
      end do
    end if

    splflo(1:nparf) = par(1:nparf)

    if ( ifs == 0) then
      call pcval(nparf,xval,tauflo,ubc,splflo)
    else if ( ifs == 1) then
      call plval(nparf,xval,tauflo,ubc,splflo)
    else if ( ifs == 2) then
      call pqval(nparf,xval,tauflo,ubc,splflo)
    else
      write ( *, * ) ' '
      write ( *, * ) 'FlowBC - Fatal error!'
      write ( *, * ) '  Illegal value of IFS = ',ifs
      stop
    end if

    vbc = 0.0D+00
!
!  Inflow points for cavity2 have the form (X,0) or (X,YRANGE).
!  NPARF must be at least 2.
!
  else if ( s_eqi ( region,'cavity2')) then

    if ( nparf == 2) then
      tauflo(1) = xrange/2.0D+00
    else
      do i = 1, nparf/2
        tauflo(i) = xrange*real ( i - 1, kind = 8 )/dble(nparf/2-1)
      end do
    end if

    if ( yval == 0.0D+00 ) then
 
      do i = 1, nparf/2
        splflo(i) = par(i)
      end do

    else if ( yval == yrange ) then
 
      do i = 1, nparf/2
        splflo(i) = par(i+nparf/2)
      end do

    end if

    if ( ifs == 0) then
      call pcval(nparf/2,xval,tauflo,ubc,splflo)
    else if ( ifs == 1) then
      call plval(nparf/2,xval,tauflo,ubc,splflo)
    else if ( ifs == 2) then
      call pqval(nparf/2,xval,tauflo,ubc,splflo)
    else
      write ( *, * ) ' '
      write ( *, * ) 'FlowBC - Fatal error!'
      write ( *, * ) '  Illegal value of IFS = ',ifs
      stop
    end if

    vbc = 0.0D+00
!
!  Inflow points for the channel have the form (0,Y).
!
!  NPARF must be at least 2, which specifies values at
!  (0,0) and (0,YRANGE).
!
  else if ( s_eqi ( region,'channel')) then

    do i = 1, nparf
      tauflo(i) = yrange*real ( i - 1, kind = 8 )/dble(nparf-1)
    end do

    splflo(1:nparf) = par(1:nparf)

    if ( ifs == 0) then
      call pcval(nparf,yval,tauflo,ubc,splflo)
    else if ( ifs == 1) then
      call plval(nparf,yval,tauflo,ubc,splflo)
    else if ( ifs == 2) then
      call pqval(nparf,yval,tauflo,ubc,splflo)
    else
      write ( *, * ) ' '
      write ( *, * ) 'FlowBC - Fatal error!'
      write ( *, * ) '  Illegal value of IFS = ',ifs
      stop
    end if

    vbc = 0.0D+00
!
!  Inflow points for the step have the coordinates (0,Y).
!
!  NPARF must be at least 2, which specifies values at
!  (0,0) and (0,YRANGE).
!
  else if ( s_eqi ( region,'step')) then

    do i = 1, nparf
      tauflo(i) = yrange*real ( i - 1, kind = 8 )/dble(nparf-1)
    end do

    splflo(1:nparf) = par(1:nparf)

    if ( ifs == 0) then
      call pcval(nparf,yval,tauflo,ubc,splflo)
    else if ( ifs == 1) then
      call plval(nparf,yval,tauflo,ubc,splflo)
    else if ( ifs == 2) then
      call pqval(nparf,yval,tauflo,ubc,splflo)
    else
      write ( *, * ) ' '
      write ( *, * ) 'FlowBC - Fatal error!'
      write ( *, * ) '  Illegal value of IFS = ',ifs
      stop
    end if

    vbc = 0.0D+00

  end if

  return
end
subroutine fpbcrb ( arb, indx, maxcofrb, maxnfl, nbcrb, ncofrb, &
  nelem, node, np, nx, ny, rb, xc, xrange, yc, yrange )

!*****************************************************************************80
!
!! FPBCRB evaluates the jacobian of the reduced boundary conditions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) ARB(MAXCOFRB,MAXCOFRB).
!    ARB contains the Jacobian or Picard matrix for the reduced
!    Navier Stokes system, stored as an NCOFRB by NCOFRB array.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!    If K = INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry GFL(K),
!    and an equation will be generated to determine its value.
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either because
!    the variable is specified in some other way, or because
!    (in the case of pressure), there is no coefficient associated
!    with that node.
!
!    Integer MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!    The local ordering of the nodes is suggested by this diagram:
!
!      Global nodes   Elements      NODE
!                                                     1  2  3  4  5  6
!      74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                     |  /   /|
!      73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                     |/   /  |
!      72  82  92     2   1-6-3
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!    The number of elements along a line in the X direction is
!    NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!    integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!    The number of elements along a line in the Y direction is
!    NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!    real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!         RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!         RB is generated by computing a finite element solution GFL,
!         which is saved for later reference as "GFLRB".
!         GFLRB is copied into the first column of RB.
!         Then, we compute the first NCOFRB derivatives of GFLRB with
!         respect to a parameter.  The first derivative
!         is stored in column 1 of RB, and so on. 
!         Now we compute the QR factorization of this matrix.
!         We intend that NEQNFL >> NCOFRB, and RB is a matrix with orthogonal
!         columns, so that:
!           Transpose(RB) * RB = Identity(1+NCOFRB)
!         If GFL is any set of finite element coefficients, the corresponding
!         set of reduced basis coefficients can be computed as:
!           GRB = Transpose(RB) * GFL
!         If GRB is a set of reduced basis coefficients, a corresponding
!         set of finite element coefficients can be computed as:
!           GFL = RB * GRB.
!         While it is the case that you can expand and then reduce,
!         and always get the same result, it is not the case that
!         when you reduce and then expand you get the same result!
!         It is true, for ANY GRB, that
!           GRB = Transpose(RB) * RB * GRB
!         which follows from Transpose(RB) * RB = Identity(1+NCOFRB).
!         However, for a general GFL, it is the case that
!           GFL  = /= RB * Transpose(RB) * GFL.
!         Only if GFL was generated from a reduced basis coefficient
!         vector will equality apply.  In other words, if GFL was generated
!         from a reduced basis coefficient:
!           GFL = RB * GRB
!         then
!           RB * Transpose(RB) * GFL = RB * Transpose(RB) * (RB * GRB)
!           = RB * GRB = GFL
!         so in this strictly limited case,
!           RB * Transpose(RB) = Identity(NEQNFL).
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  integer ( kind = 4 ) ibcrb
  integer ( kind = 4 ) icoffl
  integer ( kind = 4 ) icofrb
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) w
  real ( kind = 8 ) wurb
  real ( kind = 8 ) xbc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybc
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Zero out the BC rows of the matrix.
!
  arb(1:nbcrb,1:ncofrb) = 0.0D+00

  do ibcrb = 1, nbcrb
!
!  For the driven cavity, the boundary collocation points are evenly
!  spaced between the ends of the  upper boundary.
!
    xbc = xrange * dble(ibcrb)/dble(nbcrb+1)
    ybc = yrange

    icol = 1 + int ( xbc * dble ( nx - 1 ) / xrange )
    if ( nx - 1 < icol ) then
      icol = nx-1
    end if

    ielem = icol*(2*ny-2)-1
!
!  Evaluate the reduced solution UBCRB at (XBC,YBC).
!
    do icofrb = 1, ncofrb
      wurb = 0.0D+00
      do iq = 1, 6
        call qbf(ielem,iq,w,dwdx,dwdy,nelem,node,np,xc,xbc,yc,ybc)
        inode = node(iq,ielem)
        icoffl = indx(1,inode)
        wurb = wurb+rb(icoffl,icofrb)*w
      end do
      arb(ibcrb,icofrb) = wurb
    end do

  end do

  return
end
subroutine fpferb ( arb, area, grb, maxcofrb, maxelm, nbcrb, ncofrb, &
  nelem, nferb, phirb, reynld )

!*****************************************************************************80
!
!! FPFERB evaluates the reduced basis jacobian directly.
!
!  Discussion:
!
!    FPFERB computes the reduced basis jacobian without computing the
!    full basis jacobian first.
!
!    FPFERB is given
!
!      GRB, the reduced basis coefficients of an approximate solution,
!      PHIRB, the reduced basis functions, evaluated at the quadrature
!        points,
!      REYNLD, the current Reynolds number,
!
!    and computes
!
!      ARB, the reduced basis jacobian of the Navier Stokes equations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) ARB(MAXCOFRB,MAXCOFRB).
!    ARB contains the Jacobian or Picard matrix for the reduced
!    Navier Stokes system, stored as an NCOFRB by NCOFRB array.
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
! 
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXELM.
!    MAXELM is the maximum number of elements.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NFERB.
!    NFERB is the number of reduced basis coefficients that will
!    be determined via the finite element method.
!
!    real ( kind = 8 ) PHIRB(3,MAXCOFRB,15,MAXELM).
!    PHIRB contains the values of a finite element basis function
!    or its X or Y derivative, in a given element, at a given
!    quadrature point, for a particular reduced basis function.
!
!    For PHIRB(I,J,K,L), index J refers to the reduced basis
!    basis functions, for J = 0 to NCOFRB.
!
!    The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!      For the quadrature point I, and reduced basis function J,
!      in element L, PHIRB(I,J,K,L) represents the value of:
!
!        K = 1, WUrb, the finite element U velocity basis function;
!        K = 2, dWUrbdX, the X derivative of WUrb;
!        K = 3, dWUrbdY, the Y derivative of WUrb;
!        K = 4, WVrb, the finite element V velocity basis function;
!        K = 5, dWVrbdX, the X derivative of WVrb;
!        K = 6, dWVrbdY, the Y derivative of WVrb;
!        K = 7, Q, the finite element pressure basis function.
!        K = 8, dQrbdX, the X derivative of Qrb;
!        K = 9, dQrbdY, the Y derivative of Qrb.
!        K = 10, WU0rb, same as WUrb, with zero BC.
!        K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!        K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!        K = 13, WV0rb, same as WVrb, with zero BC.
!        K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!        K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the current value of the Reynolds number.
!    Normally, REYNLD is stored as PARA(NPARF+NPARB+1).
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
!
  real ( kind = 8 ) ar
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,maxelm)
  real ( kind = 8 ) dqjdx
  real ( kind = 8 ) dqjdy
  real ( kind = 8 ) dprbdx
  real ( kind = 8 ) dprbdy
  real ( kind = 8 ) durbdx
  real ( kind = 8 ) durbdy
  real ( kind = 8 ) dvrbdx
  real ( kind = 8 ) dvrbdy
  real ( kind = 8 ) dwu0dx
  real ( kind = 8 ) dwujdx
  real ( kind = 8 ) dwu0dy
  real ( kind = 8 ) dwujdy
  real ( kind = 8 ) dwv0dx
  real ( kind = 8 ) dwvjdx
  real ( kind = 8 ) dwv0dy
  real ( kind = 8 ) dwvjdy
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) icofrb
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) jcofrb
  logical s_eqi
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  real ( kind = 8 ) prb
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) urb
  real ( kind = 8 ) vrb
  real ( kind = 8 ) wu0
  real ( kind = 8 ) wuj
  real ( kind = 8 ) wv0
  real ( kind = 8 ) wvj
!
!  Zero out the FE rows of the matrix.
!
  do icofrb = nbcrb+1, nbcrb+nferb
    arb(icofrb,1:ncofrb) = 0.0D+00
  end do
!
!  Consider an element IELEM...
!
  do ielem = 1, nelem
!
!  ...and a quadrature point IQUAD...
!
    do iquad = 1, 3

      ar = area(iquad,ielem)
!
!  For the given reduced coefficients GRB, and basis functions
!  PHIRB, evaluate U, V, and P, and their spatial derivatives.
!
      call uvpqrb(dprbdx,dprbdy,durbdx,durbdy,dvrbdx,dvrbdy,grb, &
        ielem,iquad,maxcofrb,maxelm,ncofrb,phirb,prb,urb,vrb)
!
!  Consider FE reduced basis function ICOFRB.
!
      do icofrb = nbcrb+1, nbcrb+nferb

        wu0    = phirb(iquad,icofrb,10,ielem)
        dwu0dx = phirb(iquad,icofrb,11,ielem)
        dwu0dy = phirb(iquad,icofrb,12,ielem)
        wv0    = phirb(iquad,icofrb,13,ielem)
        dwv0dx = phirb(iquad,icofrb,14,ielem)
        dwv0dy = phirb(iquad,icofrb,15,ielem)
!
!  Take the derivative with respect to basis function JCOFRB.
!
        do jcofrb = 1, ncofrb

          wuj    = phirb(iquad,jcofrb,1,ielem)
          dwujdx = phirb(iquad,jcofrb,2,ielem)
          dwujdy = phirb(iquad,jcofrb,3,ielem)

          wvj    = phirb(iquad,jcofrb,4,ielem)
          dwvjdx = phirb(iquad,jcofrb,5,ielem)
          dwvjdy = phirb(iquad,jcofrb,6,ielem)

          dqjdx  = phirb(iquad,jcofrb,8,ielem)
          dqjdy  = phirb(iquad,jcofrb,9,ielem)
!
!  The horizontal momentum equations.
!
          arb(icofrb,jcofrb) = arb(icofrb,jcofrb)+ar* &
            (dwujdx*dwu0dx + dwujdy*dwu0dy+reynld &
            *(wuj*durbdx+urb*dwujdx+wvj*durbdy+vrb*dwujdy+dqjdx)*wu0)
!
!  The vertical momentum equations.
!
          arb(icofrb,jcofrb) = arb(icofrb,jcofrb)+ar* &
            (dwvjdx*dwv0dx + dwvjdy*dwv0dy +reynld &
            *(wuj*dvrbdx+urb*dwvjdx+wvj*dvrbdy+vrb*dwvjdy+dqjdy)*wv0)

        end do
      end do
    end do
  end do

  return
end
subroutine fpfl ( afl, area, eqn, gfl, indx, ldafl, maxelm, nelem, neqnfl, &
  nlband, node, np, npar, par, phifl )

!*****************************************************************************80
!
!! FPFL computes the jacobian of the residual function of the full solution.
!
!  Discussion:
!
!    The differentiated Navier Stokes functions have the form:
!
!
!  d U-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + reynld * (Wj*dUold/dx + Uold*dWj/dx+ Vold*dWj/dy) * Wi dx dy
!
!  d U-Eqn/d V-Coef:
!
!    Integral
!
!    reynld * Wj*dUold/dy * Wi dx dy
!
!  d U-Eqn/d P-Coef:
!
!    Integral
!
!    reynld * dQj/dx * Wi dx dy
!
!  d V-Eqn/d U-Coef:
!
!    Integral
!
!    reynld * Wj*dVold/dx * Wi dx dy
!
!  d V-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + reynld * (Uold*dWj/dx + Wj*dVold/dy + Vold*dWj/dy) * Wi dx dy
!
!  d V-Eqn/d P-Coef:
!
!    Integral
!
!    reynld * dQj/dy * Wi dx dy
!
!  d P-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * Qi dx dy
!
!  d P-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dy * Qi dx dy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) AFL(LDAFL,MAXNFL).
!    If Newton iteration is being carried out, AFL contains the
!    Jacobian matrix for the full system.
!    If Picard iteration is being carried out, AFL contains the
!    Picard matrix for the full system.
!
!    AFL is stored in LINPACK general band storage mode, with
!    logical dimensions (3*NLBAND+1, NEQNFL).
!
!    Where is the (I,J) entry of AFL actually stored?
!    AFL has actual storage for such an entry only if
!      -NLBAND <= I-J <= NLBAND.
!    In such a case, the (I,J) entry is actually stored in
!      AFL(I-J+2*NLBAND+1,J)
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    character ( len = 2 ) EQN(MAXNFL).
!    EQN records the "type" of each equation that will be generated, and
!    which is associated with an unknown. 
!
!    'U'  A horizontal momentum equation.
!    'UB' The condition U = 0 applied at a node on the bump.
!    'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U = 0 applied at a node on a fixed wall.
!    'U0' A dummy value of U = 0 should be set.
!
!    'V'  A vertical momentum equation.
!    'VB' The condition V = 0 applied at a node on the bump.
!    'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V = 0 applied at a node on a fixed wall.
!    'V0' A dummy value of V = 0 should be set.
!
!    'P'  A continuity equation.
!    'PB' The condition P = 0 applied at (XMAX,YMAX).
!    'P0' A dummy value of P = 0 should be set.
! 
!    real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
! 
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    integer ( kind = 4 ) LDAFL.
!    LDAFL is the first dimension of the matrix AFL as declared in
!    the main program.  LDAFL must be at least 3*NLBAND+1.
!
!    integer ( kind = 4 ) MAXELM.
!    MAXELM is the maximum number of elements allowed.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NLBAND.
!    NLBAND is the lower bandwidth of the matrix AFL.
!    The zero structure of AFL is assumed to be symmetric, and so
!    NLBAND is also the upper bandwidth of AFL.
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
! 
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
! 
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    real ( kind = 8 ) PHIFL(3,6,10,NELEM).
!    PHIFL contains the value of a finite element basis function, its
!    derivative, or other information, evaluated at the quadrature
!    points (which are the element midside nodes).
!
!    The meaning of the entry PHIFL(I,J,K,L) is as follows.
!    For the quadrature point I, and basis function J, in element L,
!    PHIFL(I,J,K,L) represents the value of:
!
!      K =  1, W, the finite element basis function for velocities;
!      K =  2, dWdX, the X derivative of W;
!      K =  3, dWdY, the Y derivative of W;
!      K =  4, Q, the finite element basis function for pressures;
!      K =  5, dQdX, the X derivative of Q;
!      K =  6, dQdY, the Y derivative of Q;
!      K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!      K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!      K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!      K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!    In particular, PHIFL(I,J,K,L) is the value of the quadratic
!    basis function W associated with local node J in element L,
!    evaluated at quadrature point I.
!
!    Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!    since there are only three linear basis functions.
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
!
  real ( kind = 8 ) afl(ldafl,neqnfl)
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(3,maxelm)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dqjdx
  real ( kind = 8 ) dqjdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwidx
  real ( kind = 8 ) dwidy
  real ( kind = 8 ) dwjdx
  real ( kind = 8 ) dwjdy
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iuse
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhor
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jprs
  integer ( kind = 4 ) jq
  integer ( kind = 4 ) jver
  logical s_eqi
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) qi
  real ( kind = 8 ) reynld
  real ( kind = 8 ) term
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) wi
  real ( kind = 8 ) wj
!
  reynld = par(npar)

  do i = 1, 3*nlband+1
    afl(i,1:neqnfl) = 0.0D+00
  end do
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1, nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1, 3

      ar = area(iquad,ielem)
!
!  Evaluate U, V and P at the IQUAD-th quadrature point.
!
      call uvpqfl(dpdx,dpdy,dudx,dudy,dvdx,dvdy,gfl,ielem,indx, &
        iquad,nelem,neqnfl,node,np,p,phifl,u,v)
!
!  Consider each node in the element.
!
      do iq = 1, 6

        ip = node(iq,ielem)

        wi = phifl(iquad,iq,1,ielem)
        dwidx = phifl(iquad,iq,2,ielem)
        dwidy = phifl(iquad,iq,3,ielem)
        qi = phifl(iquad,iq,4,ielem)

        ihor = indx(1,ip)
        iver = indx(2,ip)
        iprs = indx(3,ip)
!
!  Now compute the derivatives of the functions associated
!  with U, V and P, with respect to the coefficients associated
!  with basis vectors at each node of the element.
!
        do jq = 1, 6

          jp = node(jq,ielem)

          wj = phifl(iquad,jq,1,ielem)
          dwjdx = phifl(iquad,jq,2,ielem)
          dwjdy = phifl(iquad,jq,3,ielem)

          dqjdx = phifl(iquad,jq,5,ielem)
          dqjdy = phifl(iquad,jq,6,ielem)

          jhor = indx(1,jp)
          jver = indx(2,jp)
          jprs = indx(3,jp)
!
!  Contributions of the JHOR horizontal velocity to the U, V, and
!  P equations.
!
          if ( s_eqi ( eqn(ihor),'U')) then

            term = ar*(dwjdx*dwidx+dwjdy*dwidy+ &
              reynld*(wj*dudx+u*dwjdx+v*dwjdy)*wi)

            iuse = ihor-jhor+2*nlband+1
            afl(iuse,jhor) = afl(iuse,jhor)+term

          end if

          if ( s_eqi ( eqn(iver),'V')) then
            term = ar*(reynld*wj*dvdx*wi)
            iuse = iver-jhor+2*nlband+1
            afl(iuse,jhor) = afl(iuse,jhor)+term
          end if

          if ( 0 < iprs ) then
            if ( s_eqi ( eqn(iprs),'P')) then
              term = ar*dwjdx*qi
              iuse = iprs-jhor+2*nlband+1
              afl(iuse,jhor) = afl(iuse,jhor)+term
            end if
          end if
!
!  Contributions of the JVER vertical velocity variable to the
!  U, V and P equations.
!
          if ( s_eqi ( eqn(ihor),'U')) then
            term = ar*reynld*wj*dudy*wi
            iuse = ihor-jver+2*nlband+1
            afl(iuse,jver) = afl(iuse,jver)+term
          end if

          if ( s_eqi ( eqn(iver),'V')) then

            term = ar*(dwjdx*dwidx+dwjdy*dwidy &
              +reynld*(u*dwjdx+wj*dvdy+v*dwjdy)*wi)

            iuse = iver-jver+2*nlband+1
            afl(iuse,jver) = afl(iuse,jver)+term
          end if

          if ( 0 < iprs ) then
            if ( s_eqi ( eqn(iprs),'P')) then
              term = ar*dwjdy*qi
              iuse = iprs-jver+2*nlband+1
              afl(iuse,jver) = afl(iuse,jver)+term
            end if
          end if
!
!  Contributions of the JPRS pressure to the U and V equations.
!
          if ( 0 < jprs ) then

            if ( s_eqi ( eqn(ihor),'U')) then
              term = ar*reynld*dqjdx*wi
              iuse = ihor-jprs+2*nlband+1
              afl(iuse,jprs) = afl(iuse,jprs)+term
            end if

            if ( s_eqi ( eqn(iver),'V')) then
              term = ar*reynld*dqjdy*wi
              iuse = iver-jprs+2*nlband+1
              afl(iuse,jprs) = afl(iuse,jprs)+term
            end if

          end if

        end do
      end do
    end do
  end do
!
!  Set up the equations that enforce boundary conditions.
!
  do ip = 1, np

    ihor = indx(1,ip)
    iver = indx(2,ip)
    iprs = indx(3,ip)

    if ( s_eqi ( eqn(ihor),'UB').or. s_eqi ( eqn(ihor),'UI').or. &
         s_eqi ( eqn(ihor),'UW').or. s_eqi ( eqn(ihor),'U0')) then
      afl(2*nlband+1,ihor) = 1.0D+00
    end if

    if ( s_eqi ( eqn(iver),'VB').or.s_eqi ( eqn(iver),'VI').or. &
         s_eqi ( eqn(iver),'VW').or.s_eqi ( eqn(iver),'V0')) then
      afl(2*nlband+1,iver) = 1.0D+00
    end if

    if ( 0 < iprs ) then
      if ( s_eqi ( eqn(iprs),'PB')) then
        afl(2*nlband+1,iprs) = 1.0D+00
      else if ( s_eqi ( eqn(iprs),'P0')) then
        afl(2*nlband+1,iprs) = 1.0D+00
      end if
    end if

  end do

  return
end
subroutine fpirb ( afl, arb, area, eqn, gflrb, grb, indx, ldafl, &
  maxcofrb, maxelm, maxnfl, nbcrb, ncofrb, nelem, neqnfl, nferb, &
  nlband, node, np, npar, nx, ny, par, phifl, rb, xc, xrange, yc, &
  yrange )

!*****************************************************************************80
!
!! FPIRB computes the reduced basis jacobian using the indirect method.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
!
!  For some reason, I couldn't set AFL as a local array.
!
  real ( kind = 8 ) afl(ldafl,maxnfl)
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,nelem)
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gflrb(neqnfl)
!
!  FORTRAN 90 temporary array.
!
  real ( kind = 8 ) gfltmp(neqnfl)
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  logical s_eqi
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Zero out ARB.
!
  arb(1:ncofrb,1:ncofrb) = 0.0D+00
!
!  Get the partial derivative of the boundary conditions.
!
  call fpbcrb(arb,indx,maxcofrb,maxnfl,nbcrb,ncofrb, &
    nelem,node,np,nx,ny,rb,xc,xrange,yc,yrange)
!
!  Recover the equivalent full basis coefficients GFLTMP from
!  the reduced basis coefficients GRB.
!
  call grb2fl(gfltmp,gflrb,grb,maxnfl,ncofrb,neqnfl,rb)
!
!  Get the jacobian FPFL for the full coefficient set.
!
  call fpfl(afl,area,eqn,gfltmp,indx,ldafl,maxelm,nelem,neqnfl, &
    nlband,node,np,npar,par,phifl)
!
!  Zero out all rows except for U and V momentum equations.
!
  do i = 1, np
    ieqn = indx(1,i)
    if ( .not. s_eqi ( eqn(ieqn),'U')) then
      jlo = max(1,ieqn-nlband)
      jhi = min(neqnfl,ieqn+nlband)
      do j = jlo, jhi
        afl(ieqn-j+2*nlband+1,j) = 0.0D+00
      end do
    end if
    ieqn = indx(2,i)
    if ( .not. s_eqi ( eqn(ieqn),'V')) then
      jlo = max(1,ieqn-nlband)
      jhi = min(neqnfl,ieqn+nlband)
      do j = jlo, jhi
        afl(ieqn-j+2*nlband+1,j) = 0.0D+00
      end do
    end if
    ieqn = indx(3,i)
    if ( 0 < ieqn ) then
      jlo = max(1,ieqn-nlband)
      jhi = min(neqnfl,ieqn+nlband)
      do j = jlo, jhi
        afl(ieqn-j+2*nlband+1,j) = 0.0D+00
      end do
    end if
  end do
!
!  Compute the FE portion of the jacobian,
!    FPRB = QT * FPFL * Q.
!
  do i = nbcrb+1,nbcrb+nferb
    do j = 1, neqnfl
      do k = 1, neqnfl
       do l = 1, ncofrb
         if ( -nlband <= j-k.and.j-k <= nlband) then
           arb(i,l) = arb(i,l)+rb(j,i)*afl(j-k+2*nlband+1,k)*rb(k,l)
         end if
       end do
     end do
    end do
  end do

  return
end
subroutine fprb(arb,area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
  nelem,nferb,node,np,nx,ny,phirb,rb,reynld,xc,xrange,yc,yrange)

!*****************************************************************************80
!
!! FPRB evaluates the reduced basis jacobian directly.
!
!  Discussion:
!
!    Direct evaluation means that the full basis jacobian is NOT computed.
!
!    FPRB is given
!
!      PAR, the current parameter values,
!      GRB, the reduced basis coefficients of an approximate solution,
!      PHIRB, the reduced basis functions, evaluated at the quadrature points,
!
!    and computes
!
!      ARB, the reduced basis jacobian of the Navier Stokes equations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 1996.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) ARB(MAXNRB,MAXNRB).
!    ARB contains the Jacobian matrix for the reduced basis system,
!    stored as an NCOFRB by NCOFRB array.
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!  GRB    Input, real ( kind = 8 ) GRB(NCOFRB).
!         GRB contains the reduced basis coefficients of the current
!         estimate of the state solution.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB.
!         NCOFRB is the number of basis functions, reduced state equations and
!         coefficients in the reduced basis system.
!
!  PHIRB  Input, real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    Input, real ( kind = 8 ) REYNLD.
!    The current value of the Reynolds number parameter.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange

  call fpbcrb(arb,indx,maxcofrb,maxnfl,nbcrb,ncofrb, &
    nelem,node,np,nx,ny,rb,xc,xrange,yc,yrange)

  call fpferb(arb,area,grb,maxcofrb,maxelm,nbcrb,ncofrb, &
    nelem,nferb,phirb,reynld)

  return
end
subroutine fxbcrb(grb,indx,maxcofrb,maxnfl,nbcrb,ncofrb,nelem,node,np,npar, &
  nparf,nx,ny,par,rb,resrb,tauflo,xc,xrange,yc,yrange)

!*****************************************************************************80
!
!! FXBCRB evaluates the reduced boundary conditions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 1996.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
! 
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!
!    The number of elements along a line in the X direction is
!    NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!    integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!
!    The number of elements along a line in the Y direction is
!    NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!    RB is generated by computing a finite element solution GFL,
!    which is saved for later reference as "GFLRB".
!    GFLRB is copied into the first column of RB.
!    Then, we compute the first NCOFRB derivatives of GFLRB with
!    respect to a parameter.  The first derivative
!    is stored in column 1 of RB, and so on. 
!
!    real ( kind = 8 ) RESRB(NCOFRB).
!    RESRB contains the residual in the reduced basis equations,
!    for the parameter values PAR and reduced basis coefficients GRB.
!
!    real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow.  There are NPARF of them, because the end
!    values of the spline are constrained to have particular
!    values.
!
!    real ( kind = 8 ) XC(NP).
!    XC contains the X coordinates of the nodes.
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) ibcrb
  integer ( kind = 4 ) icoffl
  integer ( kind = 4 ) icofrb
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) inode
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) ubcrb
  real ( kind = 8 ) w
  real ( kind = 8 ) wurb
  real ( kind = 8 ) xbc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybc
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
  do ibcrb = 1, nbcrb
!
!  These (X,Y) values are only valid for the driven cavity.
!  You should pass REGION in here to sort it out.
!
    xbc = tauflo(ibcrb)
    ybc = yrange

    icol = 1+xbc*dble(nx-1)/xrange
    if ( nx - 1 < icol ) then
      icol = nx-1
    end if

    ielem = icol*(2*ny-2)-1
!
!  Evaluate the reduced solution UBCRB at (XBC,YBC).
!
    ubcrb = 0.0D+00
    do icofrb = 1, ncofrb
      wurb = 0.0D+00
      do iq = 1, 6
        call qbf(ielem,iq,w,dwdx,dwdy,nelem,node,np,xc,xbc,yc,ybc)
        inode = node(iq,ielem)
        icoffl = indx(1,inode)
        wurb = wurb+rb(icoffl,icofrb)*w
      end do
      ubcrb = ubcrb+grb(icofrb)*wurb
    end do
!
!  Set the function value.
!
    resrb(ibcrb) = ubcrb-par(ibcrb)

  end do

  return
end
subroutine fxferb(area,grb,maxcofrb,maxelm,nbcrb,ncofrb,nelem, &
  nferb,phirb,resrb,reynld)

!*****************************************************************************80
!
!! FXFERB evaluates the finite element portion of the reduced function.
!
!  FXFERB is given
!    GRB, the reduced basis coefficients of an approximate solution,
!    PHIRB, the reduced basis functions, evaluated at the quadrature
!      points,
!  and computes
!    RESRB, the reduced basis residual of the Navier Stokes 
!      equations.
!
!  The reduced discretized Navier Stokes equations have the form:
!
!    Integral
!
!      dUrb/dx * dWu0/dx + dUrb/dy * dWu0dy
!    + reynld * (Urb*dUrb/dx + Vrb*dUrb/dy + dPrb/dx) * Wu0 dx dy = 0
!
!    Integral
!
!      dVrb/dx * dWv0/dx + dVrb/dy * dWv0/dy
!    + reynld * (Urb*dVrb/dx + Vrb*dVrb/dy + dPrb/dy) * Wv0 dx dy = 0
!
!  Here, WU0 and WV0 are the reduced basis functions for U and V
!  assuming homogeneous boundary conditions.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!  GRB    Input, real ( kind = 8 ) GRB(NCOFRB).
!         GRB contains the reduced basis coefficients of the current
!         estimate of the state solution.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB.
!         NCOFRB is the number of basis functions, reduced state equations and
!         coefficients in the reduced basis system.
!
!  PHIRB  real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!  RESRB  Output, real ( kind = 8 ) RESRB(NCOFRB).
!         RESRB contains the residual in the reduced basis equations,
!         for the given parameter values and reduced basis coefficients GRB.
!
!    Input, real ( kind = 8 ) REYNLD.
!    The current value of the Reynolds number parameter.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) ncofrb
!
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwu0dx
  real ( kind = 8 ) dwu0dy
  real ( kind = 8 ) dwv0dx
  real ( kind = 8 ) dwv0dy
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) icofrb
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) iquad
  logical s_eqi
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  real ( kind = 8 ) p
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) wu0
  real ( kind = 8 ) wv0
!
  do icofrb = nbcrb+1, nbcrb+nferb
    resrb(icofrb) = 0.0D+00
  end do
!
!  Consider an element IELEM...
!
  do ielem = 1, nelem
!
!  ...and a quadrature point IQUAD...
!
    do iquad = 1, 3

      ar = area(iquad,ielem)
!
!  For the given reduced coefficients GRB, and basis functions
!  PHIRB, evaluate U, V, and P, and their spatial derivatives.
!
      call uvpqrb(dpdx,dpdy,dudx,dudy,dvdx,dvdy,grb, &
        ielem,iquad,maxcofrb,maxelm,ncofrb,phirb,p,u,v)
!
!  Now consider the residual associated with each finite element
!  reduced basis function.
!
      do icofrb = nbcrb+1, nbcrb+nferb

        wu0 = phirb(iquad,icofrb,10,ielem)
        dwu0dx = phirb(iquad,icofrb,11,ielem)
        dwu0dy = phirb(iquad,icofrb,12,ielem)
        wv0 = phirb(iquad,icofrb,13,ielem)
        dwv0dx = phirb(iquad,icofrb,14,ielem)
        dwv0dy = phirb(iquad,icofrb,15,ielem)

        resrb(icofrb) = resrb(icofrb)+ar* &
             (dudx*dwu0dx+dudy*dwu0dy &
             +reynld*(u*dudx+v*dudy+dpdx)*wu0 &
             +dvdx*dwv0dx+dvdy*dwv0dy &
             +reynld*(u*dvdx+v*dvdy+dpdy)*wv0)

      end do
    end do
  end do

  return
end
subroutine fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar,nparf,par, &
  phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)

!*****************************************************************************80
!
!! FXFL computes the residual of the full Navier Stokes equations.
!
!  The discretized Navier Stokes equations have the form:
!
!    Integral
!
!      dU/dx * dW/dx + dU/dy * dW/dy
!    + reynld * (U*dU/dx + V*dU/dy + dP/dx) * W dx dy = 0
!
!    Integral
!
!      dV/dx * dW/dx + dV/dy * dW/dy
!    + reynld * (U*dV/dx + V*dV/dy + dP/dy) * W dx dy = 0
!
!    Integral
!
!      (dU/dx + dV/dy) * Q dx dy = 0
!
!  Here W is a basis function for U and V, and Q is a basis
!  function for P.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 1996.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    character ( len = 2 ) EQN(MAXNFL).
!    EQN records the "type" of each equation that will be generated, and
!    which is associated with an unknown. 
!
!    'U'  A horizontal momentum equation.
!    'UB' The condition U = 0 applied at a node on the bump.
!    'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U = 0 applied at a node on a fixed wall.
!    'U0' A dummy value of U = 0 should be set.
!
!    'V'  A vertical momentum equation.
!    'VB' The condition V = 0 applied at a node on the bump.
!    'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V = 0 applied at a node on a fixed wall.
!    'V0' A dummy value of V = 0 should be set.
!
!    'P'  A continuity equation.
!    'PB' The condition P = 0 applied at (XMAX,YMAX).
!    'P0' A dummy value of P = 0 should be set.
!
!    real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
!
!    integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!    If K = INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry GFL(K),
!    and an equation will be generated to determine its value.
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either because
!    the variable is specified in some other way, or because
!    (in the case of pressure), there is no coefficient associated
!    with that node.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    The local ordering of the nodes is suggested by this diagram:
!
!      Global nodes   Elements      NODE
!                                                     1  2  3  4  5  6
!      74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                     |  /   /|
!      73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                     |/   /  |
!      72  82  92     2   1-6-3
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    real ( kind = 8 ) PHIFL(3,6,10,NELEM).
!    PHIFL contains the value of a finite element basis function, its
!    derivative, or other information, evaluated at the quadrature
!    points (which are the element midside nodes).
!
!    The meaning of the entry PHIFL(I,J,K,L) is as follows.
!    For the quadrature point I, and basis function J, in element L,
!    PHIFL(I,J,K,L) represents the value of:
!
!      K =  1, W, the finite element basis function for velocities;
!      K =  2, dWdX, the X derivative of W;
!      K =  3, dWdY, the Y derivative of W;
!      K =  4, Q, the finite element basis function for pressures;
!      K =  5, dQdX, the X derivative of Q;
!      K =  6, dQdY, the Y derivative of Q;
!      K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!      K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!      K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!      K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!    In particular, PHIFL(I,J,K,L) is the value of the quadratic
!    basis function W associated with local node J in element L,
!    evaluated at quadrature point I.
!
!    Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!    since there are only three linear basis functions.
!
!    character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottome, with tangential velocity specifications
!    there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    real ( kind = 8 ) RESFL(NEQNFL).
!    RESFL contains the residual in the full basis equations.
!
!    real ( kind = 8 ) SPLFLO(NPARF).
!    SPLFLO contains the spline coefficients for the inflow.
!
!    real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow.  There are NPARF of them, because the end
!    values of the spline are constrained to have particular
!    values.
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwidx
  real ( kind = 8 ) dwidy
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) qi
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) reynld
  logical s_eqi
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) u
  real ( kind = 8 ) ubc
  real ( kind = 8 ) v
  real ( kind = 8 ) vbc
  real ( kind = 8 ) wi
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xval
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
  real ( kind = 8 ) yval
!
  reynld = par(npar)

  resfl(1:neqnfl) = 0.0D+00
!
!  Consider an element.
!
  do ielem = 1, nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1, 3

      ar = area(iquad,ielem)
!
!  Evaluate U, V and P at the IQUAD-th quadrature point.
!
      call uvpqfl(dpdx,dpdy,dudx,dudy,dvdx,dvdy,gfl,ielem,indx, &
        iquad,nelem,neqnfl,node,np,p,phifl,u,v)
!
!  Look at nearby basis functions.
!
      do iq = 1, 6

        ip = node(iq,ielem)

        wi = phifl(iquad,iq,1,ielem)
        dwidx = phifl(iquad,iq,2,ielem)
        dwidy = phifl(iquad,iq,3,ielem)
        qi = phifl(iquad,iq,4,ielem)
!
!  The horizontal velocity equations.
!
        ihor = indx(1,ip)

        if ( s_eqi ( eqn(ihor),'U')) then

          resfl(ihor) = resfl(ihor)+ar*(dudx*dwidx + dudy*dwidy &
               +reynld*(u*dudx+v*dudy+dpdx)*wi )

        else if ( s_eqi ( eqn(ihor),'UB')) then

          resfl(ihor) = gfl(ihor)

        else if ( s_eqi ( eqn(ihor),'UI')) then

          yval = yc(ip)
          call flowbc(ifs,npar,nparf,par,region,splflo,tauflo, &
            ubc,vbc,xrange,xval,yrange,yval)

          resfl(ihor) = gfl(ihor)-ubc

        else if ( s_eqi ( eqn(ihor),'UW')) then

          resfl(ihor) = gfl(ihor)

        else if ( s_eqi ( eqn(ihor),'U0')) then

          resfl(ihor) = gfl(ihor)

        end if
!
!  The vertical velocity equations.
!
        iver = indx(2,ip)

        if ( s_eqi ( eqn(iver),'V')) then

          resfl(iver) = resfl(iver)+ar*(dvdx*dwidx + dvdy*dwidy &
               +reynld*(u*dvdx+v*dvdy+dpdy)*wi )

        else if ( s_eqi ( eqn(iver),'VB')) then

          resfl(iver) = gfl(iver)

        else if ( s_eqi ( eqn(iver),'VI')) then

          yval = yc(ip)
          call flowbc(ifs,npar,nparf,par,region,splflo,tauflo, &
              ubc,vbc,xrange,xval,yrange,yval)

          resfl(iver) = gfl(iver)-vbc

        else if ( s_eqi ( eqn(iver),'VW')) then

          resfl(iver) = gfl(iver)

        else if ( s_eqi ( eqn(iver),'V0')) then

          resfl(iver) = gfl(iver)

        end if
!
!  The pressure equations.
!
        iprs = indx(3,ip)
        if ( 0 < iprs ) then
          if ( s_eqi ( eqn(iprs),'P')) then
            resfl(iprs) = resfl(iprs)+ar*(dudx+dvdy)*qi
          else if ( s_eqi ( eqn(iprs),'PB')) then
            resfl(iprs) = gfl(iprs)
          else if ( s_eqi ( eqn(iprs),'P0')) then
            resfl(iprs) = gfl(iprs)
          end if
        end if

      end do
    end do
  end do

  return
end
subroutine fxfl2rb(grb,indx,maxcofrb,maxnfl,nbcrb,ncofrb,nelem, &
  neqnfl,nferb,node,np,npar,nparf,nx,ny,par,rb,resfl,resrb, &
  tauflo,xc,xrange,yc,yrange)

!*****************************************************************************80
!
!! FXFL2RB projects a full residual into a reduced residual.
!
!  Discussion:
!
!    The relationship used is
!
!      RESRB = Q^T * RESFL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) RESFL(NEQNFL), the function value
!    in the full system.
!
!    Output, real ( kind = 8 ) RESRB(NCOFRB), the function value
!    in the reduced system.
!
!    Input, integer ( kind = 4 ) MAXNFL, the maximum value of NEQNFL, used as
!    the leading dimension of RB.
!
!    Input, integer ( kind = 4 ) NCOFRB, the number of coefficients for the
!    reduced system.
!
!    Input, integer ( kind = 4 ) NEQNFL, the number of coefficients for the
!    full system.
!
!    Input, real ( kind = 8 ) rb(maxnfl,ncofrb).
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Compute the boundary conditions directly from GFLBC.
!
  call fxbcrb(grb,indx,maxcofrb,maxnfl,nbcrb,ncofrb,nelem, &
    node,np,npar,nparf,nx,ny,par,rb,resrb,tauflo,xc,xrange,yc,yrange)
!
!  Multiply the second half of the vector RESFL by the second half of
!  RB transpose.
!
  do i = nbcrb+1, nbcrb+nferb
    resrb(i) = 0.0D+00
    do j = 1, neqnfl
      resrb(i) = resrb(i) + rb(j,i) * resfl(j)
    end do
  end do

  return
end
subroutine fxirb ( area, eqn, gflrb, grb, ifs, indx, maxcofrb, maxnfl, &
  nbcrb, ncofrb, nelem, neqnfl, nferb, node, np, npar, nparf, nx, ny, &
  par, phifl, rb, region, resrb, splflo, tauflo, xc, xrange, yc, &
  yrange )

!***********************************************************************
!
!! FXIRB indirectly computes the reduced basis residual.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!  EQN    Input, character ( len = 2 ) EQN(NEQNFL).
!         EQN records the "type" of each equation that will be generated, and
!         which is associated with an unknown.  Note that most boundary
!         conditions do not result in an equation.  The current values are:
!
!         'U'  The horizontal momentum equation.
!         'UB' The condition U = 0 applied at a node on the bump.
!         'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!         'UW' The condition U = 0 applied at a node on a fixed wall.
!
!         'V'  The vertical momentum equation.
!         'VB' The condition V = 0 applied at a node on the bump.
!         'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!         'VW' The condition V = 0 applied at a node on a fixed wall.
!
!         'P'  The continuity equation.
!         'PB' The condition P = 0 applied at (XMAX,YMAX).
!
!  GFL    Input/output, real ( kind = 8 ) GFL(NEQNFL).
!
!         GFL must contain on input the coefficients
!         for the full basis system that are equivalent to GRB.
!
!  GFLRB  Input, real ( kind = 8 ) GFLRB(NEQNFL), the full basis coefficients
!         of the solution at which the reduced basis was generated.
!
!  GRB    Input, real ( kind = 8 ) GRB(NCOFRB).
!         The coefficients for the reduced basis system.
!
!  IFS    Input, integer ( kind = 4 ) IFS.
!         1, the inflow is modeled by C0 linear splines.
!         2, the inflow is modeled by C0 quadratic splines.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  MAXNFL Input, integer ( kind = 4 ) MAXNFL, the maximum number of equations in the
!         full system.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of equations in the full system.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB, the number of equations in the reduced
!         system.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,MAXELM).
!
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!  
!  NPAR   Input, integer ( kind = 4 ) NPAR.
!
!         The number of parameters.  NPAR = NPARF + NPARB + 1.
!
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!
!         The number of parameters associated with the position and
!         shape of the bump.
!
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1. 
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!
!         PAR is the current estimate for the parameters.
!
!  PHIFL  Input, real ( kind = 8 ) PHIFL(3,6,10,NELEM). 
!
!         PHIFL contains the value of a finite element basis function, its
!         derivative, or other information, evaluated at the quadrature
!         points.
!
!         The meaning of the entry PHIFL(I,J,K,L) is as follows. 
!         For the quadrature point I, and basis function J, in element L,
!         PHIFL(I,J,K,L) represents the value of:
!
!           K =  1, W, the finite element basis function for velocities;
!           K =  2, dWdX, the X derivative of W;
!           K =  3, dWdY, the Y derivative of W;
!           K =  4, Q, the finite element basis function for pressures;
!           K =  5, dQdX, the X derivative of Q;
!           K =  6, dQdY, the Y derivative of Q;
!           K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!           K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!           K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!           K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!         In particular, PHIFL(I,J,K,L) is the value of the quadratic
!         basis function W associated with local node J in element L,
!         evaluated at quadrature point I.
!
!         Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!         since there are only three linear basis functions.
!
!    Input, real ( kind = 8 ) RB(MAXNFL,NCOFRB), the columns of RB
!    contain the orthonormal reduced basis vectors.
!
!    Input, character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottome, with tangential velocity specifications
!    there.
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    Output, real ( kind = 8 ) RESFL(NEQNFL), the residual in the
!    full basis equations.
!
!    Output, real ( kind = 8 ) RESRB(NCOFRB), the residual in the
!    reduced basis equations, evaluated at the coefficient
!    vector GRB.
!
!    Input, real ( kind = 8 ) SPLFLO(NPARF).
!    SPLFLO contains the spline coefficients for the inflow.
!
!    Input, real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow.  There are NPARF of them, because the end
!    values of the spline are constrained to have particular
!    values.
!
!    Input, real ( kind = 8 ) XC(NP).
!    The X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YC(NP).
!    The Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) area(3,nelem)
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gflrb(neqnfl)
  real ( kind = 8 ) gfltmp(neqnfl)
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ieqn
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) indx(3,np)
  logical s_eqi
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  character ( len = 20 ) region
  real ( kind = 8 ) resfltmp(neqnfl)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Recover the equivalent full basis coefficients GFLTMP from
!  the reduced basis coefficients GRB.
!
  call grb2fl ( gfltmp, gflrb, grb, maxnfl, ncofrb, neqnfl, rb )
!
!  Evaluate the residual RESFLTMP at GFLTMP.
!
  call fxfl(area,eqn,gfltmp,ifs,indx,nelem,neqnfl,node,np,npar, &
    nparf,par,phifl,region,resfltmp,splflo,tauflo,xrange,yc,yrange)
!
!  Zero out all residuals except for U and V momentum equations.
!
  do i = 1, np
    ieqn = indx(1,i)
    if ( .not. s_eqi ( eqn(ieqn), 'U') ) then
      resfltmp(ieqn) = 0.0D+00
    end if
    ieqn = indx(2,i)
    if ( .not. s_eqi ( eqn(ieqn), 'V' ) ) then
      resfltmp(ieqn) = 0.0D+00
    end if
    ieqn = indx(3,i)
    if ( 0 < ieqn ) then
      resfltmp(ieqn) = 0.0D+00
    end if
  end do
!
!  Project the residual RESFLTMP back into the reduced space,
!  arriving at RESRB.
!
  call fxfl2rb(grb,indx,maxcofrb,maxnfl,nbcrb,ncofrb,nelem, &
    neqnfl,nferb,node,np,npar,nparf,nx,ny,par,rb,resfltmp,resrb, &
    tauflo,xc,xrange,yc,yrange)

  return
end
subroutine fxrb(area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
  nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb,resrb,reynld, &
  tauflo,xc,xrange,yc,yrange)

!*****************************************************************************80
!
!! FXRB evaluates the reduced boundary and finite element equations.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    31 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  AREA   Input, real ( kind = 8 ) AREA(3,NELEM).
!
!         AREA contains a common factor multiplying the term associated
!         with a quadrature point in a given element, namely,
!
!           AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!
!         or, if the element is isoperimetric,
!
!           AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!
!         Here Ar(IELEM) represents the area of element IELEM.
!
!  GFLBC  real ( kind = 8 ) GFLBC(NEQNFL).
!         GFLBC contains the current full solution, or, in fact,
!         ANY full solution, which satisfies the boundary conditions.
!
!  GRB    real ( kind = 8 ) GRB(NCOFRB).
!         GRB contains the reduced basis coefficients of the current
!         estimate of the state solution.
!
!  INDX   integer ( kind = 4 ) INDX(3,NP).
!         INDX(I,J) contains, for each node J, the global index of U,
!         V and P at that node, or 0 or a negative value.  The global
!         index of U, V, or P is the index of the coefficient vector
!         that contains the value of the finite element coefficient
!         associated with the corresponding basis function at the
!         given node.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  MAXCOFRB
!         Integer MAXCOFRB.
!         MAXCOFRB is the maximum legal value for NCOFRB, the number
!         of coefficients used to specify a particular reduced basis
!         solution.
!
!  MAXNFL integer MAXNFL.
!         MAXNFL is the maximum number of equations or coefficients allowed
!         for the full system.  MAXNFL must be used instead of NEQNFL as
!         the leading dimension of certain multi-dimensional arrays.
!
!  NBCRB  integer ( kind = 4 ) NBCRB.
!         NBCRB is the number of independent boundary condition
!         vectors used for the reduced basis.  NBCRB is normally
!         at least 1, and must be no more than MAXBCRB.
!
!  NCOFRB integer NCOFRB.
!         NCOFRB is the number of coefficients needed to determine
!         a particular reduced basis function.
!         NCOFRB is the sum of NBCRB and NFERB.
!
!  NELEM  integer ( kind = 4 ) NELEM.
!         NELEM is the number of elements.
!         NELEM can be determined as 2*(NX-1)*(NY-1).
!
!  NEQNFL integer NEQNFL.
!         NEQNFL is the number of equations (and coefficients) in the full
!         finite element system.
!
!  NODE   integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!         NODE(I,J) contains, for an element J, the global index of
!         the node whose local number in J is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!           Global nodes   Elements      NODE
!                                                          1  2  3  4  5  6
!           74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                          |  /   /|
!           73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                          |/   /  |
!           72  82  92     2   1-6-3
!
!  NP     integer ( kind = 4 ) NP.
!         NP is the number of nodes used to define the finite element mesh.
!         Typically, the mesh is generated as a rectangular array, with
!         an odd number of nodes in the horizontal and vertical directions.
!         The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!  NX     integer ( kind = 4 ) NX.
!         NX controls the spacing of nodes and elements in
!         the X direction.  There are 2*NX-1 nodes along various
!         lines in the X direction.
!
!         The number of elements along a line in the X direction is
!         NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!  NY     integer ( kind = 4 ) NY.
!         NY controls the spacing of nodes and elements in
!         the Y direction.  There are 2*NY-1 nodes along various
!         lines in the Y direction.
!
!         The number of elements along a line in the Y direction is
!         NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!  PHIRB  real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!  RB     real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!         RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!         RB is generated by computing a finite element solution GFL,
!         which is saved for later reference as "GFLRB".
!         GFLRB is copied into the first column of RB.
!         Then, we compute the first NCOFRB derivatives of GFLRB with
!         respect to a parameter.  The first derivative
!         is stored in column 1 of RB, and so on. 
!         Now we compute the QR factorization of this matrix.
!         We intend that NEQNFL >> NCOFRB, and RB is a matrix with orthogonal
!         columns, so that:
!           Transpose(RB) * RB = Identity(1+NCOFRB)
!         If GFL is any set of finite element coefficients, the corresponding
!         set of reduced basis coefficients can be computed as:
!           GRB = Transpose(RB) * GFL
!         If GRB is a set of reduced basis coefficients, a corresponding
!         set of finite element coefficients can be computed as:
!           GFL = RB * GRB.
!         While it is the case that you can expand and then reduce,
!         and always get the same result, it is not the case that
!         when you reduce and then expand you get the same result!
!         It is true, for ANY GRB, that
!           GRB = Transpose(RB) * RB * GRB
!         which follows from Transpose(RB) * RB = Identity(1+NCOFRB).
!         However, for a general GFL, it is the case that
!           GFL  = /= RB * Transpose(RB) * GFL.
!         Only if GFL was generated from a reduced basis coefficient
!         vector will equality apply.  In other words, if GFL was generated
!         from a reduced basis coefficient:
!           GFL = RB * GRB
!
!         then
!
!           RB * Transpose(RB) * GFL = RB * Transpose(RB) * (RB * GRB)
!           = RB * GRB = GFL
!
!         so in this strictly limited case,
!
!           RB * Transpose(RB) = Identity(NEQNFL).
!
!    real ( kind = 8 ) RESRB(NCOFRB).
!    RESRB contains the residual in the reduced basis equations,
!    for the parameter values PAR and reduced basis coefficients GRB.
!
!    Input, real ( kind = 8 ) REYNLD.
!    The current value of the Reynolds number parameter.
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
  call fxbcrb(grb,indx,maxcofrb,maxnfl,nbcrb,ncofrb,nelem, &
    node,np,npar,nparf,nx,ny,par,rb,resrb,tauflo,xc,xrange,yc,yrange)

  call fxferb(area,grb,maxcofrb,maxelm,nbcrb,ncofrb,nelem, &
    nferb,phirb,resrb,reynld)

  return
end
subroutine getgsen ( grb, gsen, icolrb, maxcofrb, nbcrb, ncofrb, &
  nsenfl, rbase )

!*****************************************************************************80
!
!! GETGSEN computes the coefficients of the sensitivity matrix S.
!
!  Discussion:
!
!    The routine uses the fact that
!      S = Q*R
!
!    Given GRB, the routine also computes the coefficients of Q. 
!
!    The calculation is simply
!      GSEN = R^(-1) * GRB
!    where R is a square upper triangular matrix.
!
!    The calculation is slightly more complicated than this, since 
!    we may have dropped some columns of Q (and hence rows and
!    columns of R.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1996
!
!  Author:
!
!    John Burkardt
! 
!  Parameters:
!
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    real ( kind = 8 ) GSEN(NBCRB+NCOFRB).
!    GSEN contains the "sensitivity coefficients".  These are simply
!    the reduced basis coefficients GRB after multiplication by
!    the inverse of RBASE, and accounting for the fact that only
!    some columns of the original set of candidate basis vectors
!    were used.
!
!    integer ( kind = 4 ) ICOLRB(MAXCOFRB).
!    ICOLRB records which columns of the initial collection of
!    candidate basis vectors were actually chosen to form the
!    reduced basis.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NSENFL.
!    NSENFL is the number of full solution sensitivities to compute,
!    counting the 0-th order sensitivity as the first one.
!
!    real ( kind = 8 ) RBASE(MAXCOFRB,MAXCOFRB).
!    RBASE is the R factor in the QR factorization of the
!    reduced basis matrix.
!    In the special case where the reduced basis matrix is
!    exactly equal to SENFL, then RBASE equals SENRB.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nsenfl
!
  real ( kind = 8 ) grb(ncofrb)
  real ( kind = 8 ) gsen(nbcrb+nsenfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) icolrb(ncofrb)
  integer ( kind = 4 ) j
  real ( kind = 8 ) rbase(maxcofrb,maxcofrb)
! 
  gsen(1:nbcrb+nsenfl) = 0.0D+00

  do i = ncofrb, 1, -1
    icol = icolrb(i)
    gsen(icol) = grb(i)
    do j = i+1, ncofrb
      gsen(icol) = gsen(icol)-rbase(i,j)*gsen(j)
    end do
    gsen(icol) = gsen(icol)/rbase(i,i)
  end do

  return
end
subroutine getbcrb ( gflrb, maxcofrb, maxnfl, nbcrb, neqnfl, rb )

!*****************************************************************************80
!
!! GETBCRB computes the vectors that will be placed into the set
!  of reduced basis vectors RB, in cases where the boundary conditions
!  depend on the parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) neqnfl
!
  real ( kind = 8 ) gflrb(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nbcrb
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
!
  do j = 1, nbcrb
!
!  Set the I-th boundary condition vector.
!  RIGHT NOW THESE ARE FAKE EQUATIONS.
!
    rb(1:neqnfl,1) = gflrb(1:neqnfl)

  end do

  return
end
subroutine getferb(icolrb,maxcofrb,maxnfl,nbcrb,ncofrb,neqnfl, &
  nferb,nsenfl,rb,rbase,senfl,senrb)

!*****************************************************************************80
!
!! GETFERB computes the finite element reduced basis vectors for RB.
!  
!  Discussion:
!
!    These vectors finish up the matrix RB.
!    The routine then orthogonalizes RB.
!
!    The routine is given:
!      NCOFRB, the number of reduced basis vectors;
!      SENFL, the full solution sensitivities of orders 1 through NCOFRB,
!    and computes
!      RB, an NEQNFL by NCOFRB orthogonal matrix, whose columns
!       were initially SENFL(0), SENFL(1), ..., SENFL(NCOFRB), and which
!       is essentially the "Q" factor of this matrix,
!      RBASE, the "R" factor in the QR factorization of the matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    17 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    integer ( kind = 4 ) ICOLRB(MAXCOFRB).
!    ICOLRB records which columns of the initial collection of
!    candidate basis vectors were actually chosen to form the
!    reduced basis.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NFERB.
!    NFERB is the number of reduced basis coefficients that will
!    be determined via the finite element method.
!
!    integer ( kind = 4 ) NSENFL.
!    NSENFL is the number of full solution sensitivities to compute,
!    counting the 0-th order sensitivity as the first one.
!
!    real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!
!    RB is generated by computing a finite element solution GFL,
!    which is saved for later reference as "GFLRB".
!    GFLRB is copied into the first column of RB.
!    Then, we compute the first NCOFRB derivatives of GFLRB with
!    respect to a parameter.  The first derivative
!    is stored in column 1 of RB, and so on. 
!
!    real ( kind = 8 ) RBASE(MAXCOFRB,MAXCOFRB).
!    RBASE is the R factor in the QR factorization of the
!    reduced basis matrix.
!
!    In the special case where the reduced basis matrix is
!    exactly equal to SENFL, then RBASE equals SENRB.
!
!    real ( kind = 8 ) SENFL(MAXNFL,MAXCOFRB).
!    Columns 1 through NSENFL of SENFL contain the sensitivities
!    of the full solution with respect to the REYNLD parameter, for
!    orders 0 through NSENFL-1.
!
!    SENFL(I,J) contains the (J-1)-th sensitivity of the I-th full unknown
!    with respect to REYNLD.
!
!    real ( kind = 8 ) SENRB(MAXCOFRB,NSENFL).
!    SENRB contains the first NSENFL order sensitivities of the
!    reduced solution with respect to the REYNLD parameter.
!
!    SENRB(I,J) contains the (J-1)-th sensitivity of the I-th reduced
!    unknown with respect to REYNLD.
!
!    SENRB is computed by premultiplying SENFL by Transpose(RB).
!      SENRB = Transpose(RB) * SENFL.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) ncofrb
!
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dtemp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icolrb(maxcofrb)
  integer ( kind = 4 ) isen
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mbcrb
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) nsenfl
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) rbase(maxcofrb,maxcofrb)
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) senrb(maxcofrb,maxcofrb)
!
!  Copy the sensitivities SENFL(:,:) into RB(:,NBCRB+1:).
!
  do i = 1, neqnfl
    do isen = 1, nsenfl
      rb(i,nbcrb+isen) = senfl(i,isen)
    end do
  end do
!
!  Initialize the R factor.
!
  do i = 1, nbcrb+nsenfl
    do j = 1, nbcrb+nsenfl
      if ( i == j) then
        rbase(i,j) = 1.0D+00
      else
        rbase(i,j) = 0.0D+00
      end if
    end do
  end do
!
!  Do a cheap sort of Gram Schmidt process to eliminate sensitivity
!  columns that are dependent on the boundary conditions or earlier
!  sensitivities.
!
  mbcrb = nbcrb
  jhi = nbcrb+nsenfl

  ncofrb = 0
  nbcrb = 0
  nferb = 0

  do j = 1, jhi
!
!  For each column of the initial RB matrix,
! 
!  ...get the Euclidean norm of the column...
!
    dtemp1 = dnrm2(neqnfl,rb(1,j),1)
!
!  ...and then subtract off the projections onto the
!  already accepted columns...
!
    do i = 1, ncofrb
      dtemp = ddot(neqnfl,rb(1,i),1,rb(1,j),1)
      rbase(i,ncofrb+1) = dtemp
      call daxpy(neqnfl,-dtemp,rb(1,i),1,rb(1,j),1)
    end do
!
!  ...then get the Euclidean norm of what is left,
!  save it in RBASE, and normalize the column.
!
    dtemp = dnrm2(neqnfl,rb(1,j),1)
    rbase(ncofrb+1,ncofrb+1) = dtemp

    do i = 1, neqnfl
      rb(i,ncofrb+1) = rb(i,j)/dtemp
    end do
!
!  Now decide whether to accept this column.  
!
    if ( dtemp /= 0.0D+00 .and. 0.00001D+00 * dtemp1 < dtemp ) then

      ncofrb = ncofrb+1

      if ( j <= mbcrb) then
        nbcrb = nbcrb+1
      else
        nferb = nferb+1
      end if
             
      icolrb(ncofrb) = j

    end if

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GETRB - Information:'
  write ( *, '(a,i6)' ) '  # of BC vectors, NBCRB  = ',nbcrb
  write ( *, '(a,i6)' ) '  # of FE vectors, NFERB  = ',nferb
  write ( *, '(a,i6)' ) '  # of RB coeffs,  NCOFRB = ',ncofrb
!
!  Compute SENRB as Transpose(Q) * SENFL.
!  NSENFL may, or may not, equal NCOFRB.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GETRB - Note:'
  write ( *, '(a)' ) '  Automatically computing SENRB.'

  do i = 1, ncofrb
    do j = 1, nsenfl
      senrb(i,j) = 0.0D+00
      do k = 1, ncofrb
        senrb(i,j) = senrb(i,j)+rb(k,i)*senfl(k,j)
      end do
    end do
  end do

  return
end
subroutine getsenfl(afl,area,eqn,gfl,indx,ipivfl,ldafl,maxcofrb,maxnfl, &
  nelem,neqnfl,nlband,node,np,npar,nsenfl,par,phifl,resfl,senfl)

!*****************************************************************************80
!
!! GETSENFL computes the matrix SENFL of sensitivity vectors.
!
!  Discussion:
!
!    The routine first saves a copy of GFL, the current solution.
!
!    Then it constructs a matrix SENFL, whose first column is d GFL/d REYNLD,
!    which must be computed by solving a linear system. 
!
!    Similarly, higher derivatives of GFL with respect to REYNLD are computed
!    and stored in successive columns of SENFL. 
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) AFL(LDAFL,MAXNFL).
!    If Newton iteration is being carried out, AFL contains the
!    Jacobian matrix for the full system.
!    If Picard iteration is being carried out, AFL contains the
!    Picard matrix for the full system.
!
!    AFL is stored in LINPACK general band storage mode, with
!    logical dimensions (3*NLBAND+1, NEQNFL).
!
!    Where is the (I,J) entry of AFL actually stored?
!    AFL has actual storage for such an entry only if
!      -NLBAND <= I-J <= NLBAND.
!    In such a case, the (I,J) entry is actually stored in
!      AFL(I-J+2*NLBAND+1,J)
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!
!    or, if the element is isoperimetric,
!
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    character ( len = 2 ) EQN(MAXNFL).
!    EQN records the "type" of each equation that will be generated, and
!    which is associated with an unknown. 
!
!    'U'  A horizontal momentum equation.
!    'UB' The condition U = 0 applied at a node on the bump.
!    'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U = 0 applied at a node on a fixed wall.
!    'U0' A dummy value of U = 0 should be set.
!
!    'V'  A vertical momentum equation.
!    'VB' The condition V = 0 applied at a node on the bump.
!    'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V = 0 applied at a node on a fixed wall.
!    'V0' A dummy value of V = 0 should be set.
!
!    'P'  A continuity equation.
!    'PB' The condition P = 0 applied at (XMAX,YMAX).
!    'P0' A dummy value of P = 0 should be set.
!
!    real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    Workspace, integer IPIVFL(NEQNFL).
!    IPIVFL is a pivot vector for the solution of the full
!    linear system.
!
!    integer ( kind = 4 ) LDAFL.
!    LDAFL is the first dimension of the matrix AFL as declared in
!    the main program.  LDAFL must be at least 3*NLBAND+1.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NLBAND.
!    NLBAND is the lower bandwidth of the matrix AFL.
!    The zero structure of AFL is assumed to be symmetric, and so
!    NLBAND is also the upper bandwidth of AFL.
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NSENFL.
!    NSENFL is the number of full solution sensitivities to compute,
!    counting the 0-th order sensitivity as the first one.
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    real ( kind = 8 ) PHIFL(3,6,10,NELEM).
!    PHIFL contains the value of a finite element basis function, its
!    derivative, or other information, evaluated at the quadrature
!    points (which are the element midside nodes).
!
!    The meaning of the entry PHIFL(I,J,K,L) is as follows.
!    For the quadrature point I, and basis function J, in element L,
!    PHIFL(I,J,K,L) represents the value of:
!
!      K =  1, W, the finite element basis function for velocities;
!      K =  2, dWdX, the X derivative of W;
!      K =  3, dWdY, the Y derivative of W;
!      K =  4, Q, the finite element basis function for pressures;
!      K =  5, dQdX, the X derivative of Q;
!      K =  6, dQdY, the Y derivative of Q;
!      K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!      K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!      K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!      K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!    In particular, PHIFL(I,J,K,L) is the value of the quadratic
!    basis function W associated with local node J in element L,
!    evaluated at quadrature point I.
!
!    Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!    since there are only three linear basis functions.
!
!    real ( kind = 8 ) RESFL(NEQNFL).
!    RESFL contains the residual in the full basis equations.
!
!    real ( kind = 8 ) SENFL(MAXNFL,MAXCOFRB).
!    Columns 1 through NSENFL of SENFL contain the sensitivities
!    of the full solution with respect to the REYNLD parameter, for
!    orders 0 through NSENFL-1.
!
!    SENFL(I,J) contains the (J-1)-th sensitivity of the I-th full unknown
!    with respect to REYNLD.
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
!
  real ( kind = 8 ) afl(ldafl,maxnfl)
  real ( kind = 8 ) area(3,nelem)
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ipivfl(maxnfl)
  integer ( kind = 4 ) isen
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nsenfl
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rpnrm
  real ( kind = 8 ) ruvnrm
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) spnrm
  real ( kind = 8 ) suvnrm
!
  reynld = par(npar)
!
!  Compute, one at a time, the columns of the RB matrix.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GETSENFL - Information:'
  write ( *, * ) '  Number of sensitivities requested, NSENFL = ',nsenfl
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    'Order     MxNorm(UVRHS) MxNorm(PRHS)  MxNorm(UVSen) MxNorm(PSen)'

  senfl(1:neqnfl,1) = gfl(1:neqnfl)

  isen = 0
  call uvpnrm(senfl(1,1),indx,neqnfl,np,spnrm,suvnrm)
  write(*,'(1x,i6,28x,2g14.6)')isen,suvnrm,spnrm

  do isen = 1, nsenfl-1
!
!  Given the current solution and lower order sensitivities
!  in SENFL, compute in RESFL the right hand side for sensitivity
!  of order ISEN.
!
    call reysen(area,eqn,indx,isen,maxcofrb,maxnfl,nelem, &
      neqnfl,node,np,phifl,resfl,reynld,senfl)
!
!  Compute the norm of this right hand side.
!
    call uvpnrm(resfl,indx,neqnfl,np,rpnrm,ruvnrm)
!
!  Solve the linear system AFL * SENFL(ISEN) = RESFL
!
    call dsolfl(afl,ldafl,neqnfl,nlband,nlband,ipivfl,resfl)
!
!  Get the norm of this new sensitivity.
!
    call uvpnrm(resfl,indx,neqnfl,np,spnrm,suvnrm)

    write(*,'(1x,i6,4g14.6)')isen,ruvnrm,rpnrm,suvnrm,spnrm
!
!  Copy the new sensitivity into the SENFL array.
!
    senfl(1:neqnfl,isen+1) = resfl(1:neqnfl)

  end do

  return
end
subroutine getsenrb(maxcofrb,maxnfl,ncofrb,neqnfl,nsenfl,rb,senfl,senrb)

!*****************************************************************************80
! 
!! GETSENRB determines the reduced sensitivities from the full sensitivities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Integer MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis solution.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NCOFRB, the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NEQNFL, the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NSENFL, the number of full solution sensitivities to compute,
!    counting the 0-th order sensitivity as the first one.
!
!  RB     real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!         RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!         RB is generated by computing a finite element solution GFL,
!         which is saved for later reference as "GFLRB".
!         GFLRB is copied into the first column of RB.
!         Then, we compute the first NCOFRB derivatives of GFLRB with
!         respect to a parameter.  The first derivative
!         is stored in column 1 of RB, and so on. 
!         Now we compute the QR factorization of this matrix.
!         We intend that NEQNFL >> NCOFRB, and RB is a matrix with orthogonal
!         columns, so that:
!           Transpose(RB) * RB = Identity(NCOFRB)
!         If GFL is any set of finite element coefficients, the corresponding
!         set of reduced basis coefficients can be computed as:
!           GRB = Transpose(RB) * GFL
!         If GRB is a set of reduced basis coefficients, a corresponding
!         set of finite element coefficients can be computed as:
!           GFL = RB * GRB.
!         While it is the case that you can expand and then reduce,
!         and always get the same result, it is not the case that
!         when you reduce and then expand you get the same result!
!         It is true, for ANY GRB, that
!           GRB = Transpose(RB) * RB * GRB
!         which follows from Transpose(RB) * RB = Identity(NCOFRB).
!         However, for a general GFL, it is the case that
!           GFL  = /= RB * Transpose(RB) * GFL.
!         Only if GFL was generated from a reduced basis coefficient
!         vector will equality apply.  In other words, if GFL was generated
!         from a reduced basis coefficient:
!           GFL = RB * GRB
!         then
!           RB * Transpose(RB) * GFL = RB * Transpose(RB) * (RB * GRB)
!           = RB * GRB = GFL
!         so in this strictly limited case,
!           RB * Transpose(RB) = Identity(NEQNFL).
!
!    real ( kind = 8 ) SENFL(MAXNFL,MAXCOFRB).
!    Columns 1 through NSENFL of SENFL contain the sensitivities
!    of the full solution with respect to the REYNLD parameter, for
!    orders 0 through NSENFL-1.
!    SENFL(I,J) contains the (J-1)-th sensitivity of the I-th full unknown
!    with respect to REYNLD.
!
!    real ( kind = 8 ) SENRB(MAXCOFRB,NSENFL).
!    SENRB contains the first NSENFL order sensitivities of the
!    reduced solution with respect to the REYNLD parameter.
!    SENRB(I,J) contains the (J-1)-th sensitivity of the I-th reduced
!    unknown with respect to REYNLD.
!    SENRB is computed by premultiplying SENFL by Transpose(RB).
!      SENRB = Transpose(RB) * SENFL.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nsenfl
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) senrb(maxcofrb,maxcofrb)
!
!  Multiply SENRB = QT * SENFL
!
  do i = 1, ncofrb
    do j = 1, nsenfl
      senrb(i,j) = 0.0D+00
      do k = 1, neqnfl
        senrb(i,j) = senrb(i,j) + rb(k,i) * senfl(k,j)
      end do
    end do
  end do

  return
end
subroutine gfl2rb ( gfl, gflrb, grb, maxnfl, ncofrb, neqnfl, rb )

!*****************************************************************************80
!
!! GFL2RB projects a full solution GFL into the reduced solution GRB.
!
!  Discussion:
!
!    The relationship used is
!
!      GRB = Q^T * (GFL-GFLRB).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 September 1996.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
!
!    Input, real ( kind = 8 ) GFLRB(NEQNFL).
!    GFLRB is the solution value at which the reduced basis was computed.
!    The corresponding parameters are PARRB.
!
!    Output, real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    Input, integer ( kind = 4 ) MAXNFL, the maximum value of NEQNFL, used as
!    the leading dimension of RB.
!
!    Input, integer ( kind = 4 ) NCOFRB, the number of coefficients for the
!    reduced system.
!
!    Input, integer ( kind = 4 ) NEQNFL, the number of coefficients for the
!    full system.
!
!    Input, real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!
  implicit none
!
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) ncofrb
!
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gflrb(neqnfl)
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) rb(maxnfl,ncofrb)
!
!  Multiply (GFL-GFLRB) by RB transpose.
!
  do i = 1, ncofrb
    grb(i) = 0.0D+00
    do j = 1, neqnfl
      grb(i) = grb(i) + rb(j,i) * ( gfl(j) - gflrb(j) )
    end do
  end do

  return
end
subroutine grb2fl ( gfl, gflrb, grb, maxnfl, ncofrb, neqnfl, rb )

!*****************************************************************************80
!
!! GRB2FL determines the full solution represented by a reduced solution.
!
!  Discussion:
!
!    The relationship used is:
!
!      GFL = GFLRB + Q*GRB.
!
!    GRB2FL is given:
!
!      NCOFRB, the number of reduced basis vectors and coefficients;
!      NEQNFL, the number of full basis vectors and coefficients;
!      GRB, the reduced basis coefficients;
!      RB, the matrix of reduced basis vectors.
!
!    and computes
!
!      GFL, the corresponding full solution.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 September 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
!
!    Input, real ( kind = 8 ) GFLRB(NEQNFL).
!    GFLRB is the solution value at which the reduced basis was computed.
!    The corresponding parameters are PARRB.
!
!    Input, real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    Input, integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    Input, integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    Input, integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!
  implicit none
!
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) neqnfl
!
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gflrb(neqnfl)
  real ( kind = 8 ) grb(ncofrb)
  real ( kind = 8 ) rb(maxnfl,ncofrb)
!
  gfl(1:neqnfl) = gflrb(1:neqnfl) &
    + matmul ( rb(1:neqnfl,1:ncofrb), grb(1:ncofrb) )

  return
end
subroutine hello ( maxnx, maxny )

!*****************************************************************************80
!
!! HELLO prints an introductory message about the program.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) MAXNX.
!    MAXNX is the maximum size of NX that the program can handle.
!
!    Input, integer ( kind = 4 ) MAXNY.
!    MAXNY is the maximum size of NY that the program can handle.
!
  implicit none
!
  integer ( kind = 4 ) maxnx
  integer ( kind = 4 ) maxny
!
!  Say hello.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARBY4'
  write ( *, '(a)' ) '  A reduced basis flow analysis code.'
  write ( *, '(a)' ) '  Last modified on 04 December 2000.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The maximum problem size is'
  write ( *, '(a,i6)' ) '    MAXNX = ', maxnx
  write ( *, '(a,i6)' ) '    MAXNY = ', maxny

  return
end
subroutine help

!*****************************************************************************80
!
!! HELP lists the interactive commands.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HELP'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Compare    Compare GFL to GFLSAV;'
  write ( *, '(a)' ) 'CostFL     Evaluate cost of GFL;'
  write ( *, '(a)' ) 'CostRB     Evaluate cost of GRB;'
  write ( *, '(a)' ) 'DetFpFL    Determinant of full jacobian;'
  write ( *, '(a)' ) 'DetFpRB    Determinant of reduced jacobian;'
  write ( *, '(a)' ) 'DifFPRB    FD estimate of reduced jacobian;'
  write ( *, '(a)' ) 'DifSenFL   FD estimate of full sensitivities;'
  write ( *, * ) 'DifSenRB   FD estimate of reduced sensitivities;'
  write ( *, * ) 'DisFil =     Name the DISPLAY output file;'
  write ( *, * ) 'DisPlot    Make DISPLAY plot file of current data;'
  write ( *, * ) 'DREY =       Set REYNLD Taylor increment;'
  write ( *, * ) 'Echo       Echo user commands;'
  write ( *, * ) 'EPSDIF =     Set finite difference increment;'
  write ( *, * ) 'Expand GRB Compute GFL = RB*GRB;'
  write ( *, * ) 'FPFL       Evaluate full jacobian;'
  write ( *, * ) 'FPIRB      Evaluate reduced jacobian indirectly;'
  write ( *, * ) 'FPRB       Evaluate reduced jacobian;'
  write ( *, * ) 'FPRB = 0     Zero out reduced jacobian;'
  write ( *, * ) 'FxFl       Evaluate full residual, FXFL(GFL);'
  write ( *, * ) 'FxIRB      Evaluate FXRB = RB^T*FX(RB*GRB), indirectly;'
  write ( *, * ) 'FxRB       Evaluate FXRB = FXrb(GRB) directly;'
  write ( *, * ) 'FxRB = 0     Set vector FXRB=0;'
  write ( *, * ) 'GetGSEN    Compute sensitivity coefficients;'
  write ( *, * ) 'GetRB      Compute reduced basis;'
  write ( *, * ) 'GetSenFL   Compute full sensitivities;'
  write ( *, * ) 'GetSenRB   Compute reduced sensitivities;'
  write ( *, * ) 'GFL =        Set current full solution,'
  write ( *, * ) '           Legal values: 0, GFLSAV, GFLTAY, TAYLOR;'
  write ( *, * ) 'GFLSAV = GFL Save current GFL value;'
  write ( *, * ) 'GFLTAY =     Set Taylor base solution,'
  write ( *, * ) '           Legal values: 0, GFL, GFLSAV;'
  write ( *, * ) 'GFLTMP =     Set temporary base solution,'
  write ( *, * ) '           Legal values: 0, GFL, GFL-GFLSAV,'
  write ( *, * ) '           GFL-GFLTAR, GFLSAV, GFLSAV-GFLTAY;'
  write ( *, * ) 'GRB(*) = *   Set an entry of GRB to a value;'
  write ( *, * ) 'GRB =        Set current reduced solution GRB,'
  write ( *, * ) '           Legal values: 0, GRBSAV, GRBTAY, TAYLOR;'
  write ( *, * ) 'GRB = (*,*,...,*)  Set individual entries of GRB;'
  write ( *, * ) 'GRBSAV =     Save a GRB value,'
  write ( *, * ) '           Legal values: 0, GRB;'
  write ( *, * ) 'GRBTAY =     Set Taylor base solution,'
  write ( *, * ) '           Legal values: 0, GRB, GRBSAV;'
  write ( *, * ) 'GridX =      Uniform, Cos, or SqrtSin;'
  write ( *, * ) 'GridY =      Uniform, Cos, or SqrtSin;'
  write ( *, * ) 'Hello      Print program version and other info;'
  write ( *, * ) 'Help       Print list of commands;'
  write ( *, * ) 'IBS =        Set bump shape option;'
  write ( *, * ) 'IBUMP =      Set bump option;'
  write ( *, * ) 'IFS =        Set inflow shape option;'
  write ( *, * ) 'IHI =        Maximum row for printout,'
  write ( *, * ) '           NCOFRB, NEQNFL, NP are legal;'
  write ( *, * ) 'IJAC =       Set Jacobian option;'
  write ( *, * ) 'ILO =        Minimum row for printout;'
  write ( *, * ) 'Init       Initialize variables;'
  write ( *, * ) 'IOPT(*) =    Specify free or fixed variables;'
  write ( *, * ) 'IWRITE =     Set level of output;'
  write ( *, * ) 'JHI =        Maximum column for printout,'
  write ( *, * ) '           NCOFRB, NEQNFL, NSENFL are legal;'
  write ( *, * ) 'JLO =        Minimum column for printout,'
  write ( *, * ) 'L2NORM *   Compute big L2 norm of *,'
  write ( *, * ) '           * = GFL, GFLSAV, GFLTAR, GFLTAY, GFLTMP;'
  write ( *, * ) 'MAXNEW =     Set number of Newton steps;'
  write ( *, * ) 'MAXOPT =     Set number of optimization steps;'
  write ( *, * ) 'MAXSIM =     Set number of simple steps;'
  write ( *, * ) 'NBCRB =      Set number of boundary conditions;'
  write ( *, * ) '           (0 or 1, right now);'
  write ( *, * ) 'NewtFL     Newton''s method applied to GFL;'
  write ( *, * ) 'NewtRB     Newton''s method applied to GRB;'
  write ( *, * ) 'NPARB =      Set number of bump parameters;'
  write ( *, * ) 'NPARF =      Set number of inflow parameters;'
  write ( *, * ) 'NSENFL =     Set number of full sensitivities;'
  write ( *, * ) 'NTAY =       Set number of Taylor vectors to use,'
  write ( *, * ) '           NTAY = NCOFRB is legal, too;'
  write ( *, * ) 'NX =         Set number of X nodes;'
  write ( *, * ) 'NY =         Set number of Y nodes;'
  write ( *, * ) 'OptDifFl   Optimize the full system, using'
  write ( *, * ) '           FD estimates for gradients;'
  write ( *, * ) 'PAR(*) =     Set a parameter;'
  write ( *, * ) 'PARTAR(*) =  Set a target parameter;'
  write ( *, * ) 'PicFL      Picard''s method applied to GFL;'
  write ( *, * ) 'PicRB      Picard''s method applied to GRB;'
  write ( *, * ) 'PrDat      Print the variable values;'
  write ( *, * ) 'PrElem     Print element data, for elements'
  write ( *, * ) '           ILO to IHI;'
  write ( *, * ) 'PrFPFL     Print full jacobian,'
  write ( *, * ) '           Equations ILO to IHI,'
  write ( *, * ) '           Variables JLO to JHI;'
  write ( *, * ) 'PrFPRB     Print reduced jacobian,'
  write ( *, * ) '           Equations ILO to IHI,'
  write ( *, * ) '           Variables JLO to JHI;'
  write ( *, * ) 'PrFXFL     Print FXFL(GFL),'
  write ( *, * ) '           nodes ILO to IHI;'
  write ( *, * ) 'PrFXFLNrm  Print norm of FXFL(GFL);'
  write ( *, * ) 'PrFXRB     Print FXRB(GRB),'
  write ( *, * ) '           equations ILO to IHI;'
  write ( *, * ) 'PrGFL      Print full solution GFL,'
  write ( *, * ) '           nodes ILO to IHI;'
  write ( *, * ) 'PrGFLNrm   Print GFL and FX(GFL) norms;'
  write ( *, * ) 'PrGRB      Print reduced solution GRB;'
  write ( *, * ) 'PrGSEN     Print sensitivity coefficients;'
  write ( *, * ) 'PrINDX     Print node/equation table INDX,'
  write ( *, * ) '           nodes ILO to IHI.'
  write ( *, * ) 'PrPar      Print current parameters;'
  write ( *, * ) 'PrParTar   Print target parameters;'
  write ( *, * ) 'PrRBase    Print the reduced basis R factor;'
  write ( *, * ) 'PrRB       Print reduced basis matrix RB,'
  write ( *, * ) '           nodes ILO to IHI,'
  write ( *, * ) '           columns JLO to JHI;'
  write ( *, * ) 'PrSenFL    Print full sensitivity matrix SENFL,'
  write ( *, * ) '           nodes ILO to IHI,'
  write ( *, * ) '           sensitivities JLO to JHI;'
  write ( *, * ) 'PrSenNrm   Print full sensitivity norms;'
  write ( *, * ) 'PrSenRB    Print reduced sensitivity matrix SENRB,'
  write ( *, * ) '           rows ILO to IHI,'
  write ( *, * ) '           columns JLO to JHI;'
  write ( *, * ) 'PrUVPGFL   Print full solution at nodes in'
  write ( *, * ) '           XMIN, YMIN, XMAX, YMAX;'
  write ( *, * ) 'PrUVPRB    Print reduced basis vectors at nodes in'
  write ( *, * ) '           XMIN, YMIN, XMAX, YMAX;'
  write ( *, * ) 'PrUVPSENFL Print sensitivity vectors at nodes in'
  write ( *, * ) '           XMIN, YMIN, XMAX, YMAX;'
  write ( *, * ) 'PrUVPGRB   Print reduced solution at nodes in'
  write ( *, * ) '           XMIN, YMIN, XMAX, YMAX;'
  write ( *, * ) 'PrXY       Print X, Y nodal coordinates,'
  write ( *, * ) '           nodes ILO to IHI;'
  write ( *, * ) 'Reduce GFL Compute GRB = RB^T * GFL;'
  write ( *, * ) 'REGION =     Cavity, Cavity2, Channel, or Step;'
  write ( *, * ) 'REYNLD =     Set REYNLD parameter;'
  write ( *, * ) 'REYTAY =     Set REYNLD parameter for Taylor;'
  write ( *, * ) '           ("REYTAY = REYNLD" is legal.)'
  write ( *, * ) 'SetGeo     Set problem geometric data;'
  write ( *, * ) 'SetLog     Set problem logical data;'
  write ( *, * ) 'Stop       Stop the program;'
  write ( *, * ) 'Target     Save current GFL as GFLTAR;'
  write ( *, * ) 'Test2      Compare full and reduced U,V,P'
  write ( *, * ) '           in elements ILO through IHI;'
  write ( *, * ) 'Test3      Compare RB*RFACT and SENFL;'
  write ( *, * ) 'Test4      Compare regular and FD full sens;'
  write ( *, * ) 'Test5      Compare RB*RFACT and old RB;'
  write ( *, * ) 'Time       Print elapsed time;'
  write ( *, * ) 'TOLNEW =     Set Newton tolerance;'
  write ( *, * ) 'TOLOPT =     Set optimization tolerance;'
  write ( *, * ) 'TOLSIM =     Set Picard tolerance;'
  write ( *, * ) 'TecFil =     Name the TECPLOT output file;'
  write ( *, * ) 'TecPlot    Make TECPLOT plot file of current data;'
  write ( *, * ) 'WATEB =      Set bump weight in cost;'
  write ( *, * ) 'WATEP =      Set pressure weight in cost;'
  write ( *, * ) 'WATEU =      Set H-velocity weight in cost;'
  write ( *, * ) 'WATEV =      Set V-velocity weight in cost;'
  write ( *, * ) 'XBL =        Set left bump X coordinate;'
  write ( *, * ) 'XBR =        Set right bump X coordinate;'
  write ( *, * ) 'XMAX =       Specify XMAX.'
  write ( *, * ) 'XMIN =       Specify XMIN.'
  write ( *, * ) 'XPROF =      Set profile X coordinate;'
  write ( *, * ) 'XRANGE =     Set width of region;'
  write ( *, * ) 'YBL =        Set left bump Y coordinate;'
  write ( *, * ) 'YBR =        Set right bump Y coordinate;'
  write ( *, * ) 'YMAX =       Specify YMAX.'
  write ( *, * ) 'YMIN =       Specify YMIN.'
  write ( *, * ) 'YRANGE =     Set height of region;'

  return
end
subroutine init(afl,arb,area,command,cost,costb,costp,costu,costv,difcof, &
  disfil,drey,epsdif,eqn,etaq,gfl,gflafl,gflrb,gflsav,gflsen,gfltar,gfltay, &
  grb,grbarb,grbsav,grbsen,grbtay,gridx,gridy,hx,hy,ibs,ibump,icolrb,ierror, &
  ifs,ihi,ijac,ilo,indx,iopt,ipar,ipivfl,ipivrb,isotri,iwrite,jhi,jlo,ldafl, &
  maxcofrb,maxelm,maxnew,maxnfl,maxnp,maxny,maxopt,maxpar,maxparb,maxparf, &
  maxsim,nbcrb,ncofrb,nelem,neqnfl,nferb,nlband,node,nodelm,np,npar,nparb, &
  nparf,npe,nprof,nsenfl,ntay,numnew,numopt,numsim,nx,ny,par,parafl,pararb, &
  pardif,parrb,parsav,parsen,partar,phifl,phirb,rb,rbase,region,resfl, &
  resflsav,resrb,reynld,reytay,senfl,senrb,splbmp,splflo,taubmp,tauflo, &
  tecfil,tolnew,tolopt,tolsim,value,wateb,watep,wateu,watev,wquad,xbl,xbr, &
  xc,xmax,xmin,xprof,xquad,xrange,xsiq,ybl,ybr,yc,ymax,ymin,yquad,yrange)

!*****************************************************************************80
!
!! INIT sets problem data to default values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) maxnp
  integer ( kind = 4 ) maxny
  integer ( kind = 4 ) maxpar
  integer ( kind = 4 ) maxparb
  integer ( kind = 4 ) maxparf
!
  real ( kind = 8 ) afl(ldafl,maxnfl)
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,maxelm)
  character ( len = 80 ) command
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) difcof(maxcofrb)
  character ( len = 30 ) disfil
  real ( kind = 8 ) drey
  real ( kind = 8 ) epsdif
  character ( len = 2 ) eqn(maxnfl)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) gfl(maxnfl)
  real ( kind = 8 ) gflafl(maxnfl)
  real ( kind = 8 ) gflrb(maxnfl)
  real ( kind = 8 ) gflsav(maxnfl)
  real ( kind = 8 ) gflsen(maxnfl)
  real ( kind = 8 ) gfltar(maxnfl)
  real ( kind = 8 ) gfltay(maxnfl)
  real ( kind = 8 ) grb(maxcofrb)
  real ( kind = 8 ) grbarb(maxcofrb)
  real ( kind = 8 ) grbsav(maxcofrb)
  real ( kind = 8 ) grbsen(maxcofrb)
  real ( kind = 8 ) grbtay(maxcofrb)
  character ( len = 20 ) gridx
  character ( len = 20 ) gridy
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) icolrb(maxcofrb)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx(3,maxnp)
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivfl(maxnfl)
  integer ( kind = 4 ) ipivrb(maxcofrb)
  integer ( kind = 4 ) isotri(maxelm)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) maxsim
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,maxelm)
  integer ( kind = 4 ) nodelm(maxnp)
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nprof(2*maxny-1)
  integer ( kind = 4 ) nsenfl
  integer ( kind = 4 ) ntay
  integer ( kind = 4 ) numnew
  integer ( kind = 4 ) numopt
  integer ( kind = 4 ) numsim
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(maxpar)
  real ( kind = 8 ) parafl(maxpar)
  real ( kind = 8 ) pararb(maxpar)
  real ( kind = 8 ) pardif ( maxpar)
  real ( kind = 8 ) parrb(maxpar)
  real ( kind = 8 ) parsav(maxpar)
  real ( kind = 8 ) parsen(maxpar)
  real ( kind = 8 ) partar(maxpar)
  real ( kind = 8 ) phifl(3,6,10,maxelm)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) rbase(maxcofrb,maxcofrb)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(maxnfl)
  real ( kind = 8 ) resflsav(maxnfl)
  real ( kind = 8 ) resrb(maxcofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) reytay
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) senrb(maxcofrb,maxcofrb)
  real ( kind = 8 ) splbmp(maxparb+2)
  real ( kind = 8 ) splflo(maxparf)
  real ( kind = 8 ) taubmp(maxparb+2)
  real ( kind = 8 ) tauflo(maxparf)
  character ( len = 30 ) tecfil
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) value
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xquad(3,maxelm)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yquad(3,maxelm)
  real ( kind = 8 ) yrange
!
!  Zero out the variables.
!
  afl(1:ldafl,1:maxnfl) = 0.0D+00
  arb(1:maxcofrb,1:maxcofrb) = 0.0D+00
  area(1:3,1:maxelm) = 0.0D+00
  command = ' '
  cost = 0.0D+00
  costb = 0.0D+00
  costp = 0.0D+00
  costu = 0.0D+00
  costv = 0.0D+00
  difcof(1:maxcofrb) = 0.0D+00
  disfil = 'display.dat'
  drey = 0.01D+00
  epsdif = 0.000001D+00
  eqn(1:maxnfl) = ' '
  etaq(1:3) = 0.0D+00
  gfl(1:maxnfl) = 0.0D+00
  gflafl(1:maxnfl) = 0.0D+00
  gflrb(1:maxnfl) = 0.0D+00
  gflsav(1:maxnfl) = 0.0D+00
  gflsen(1:maxnfl) = 0.0D+00
  gfltar(1:maxnfl) = 0.0D+00
  gfltay(1:maxnfl) = 0.0D+00
  grb(1:maxcofrb) = 0.0D+00
  grbarb(1:maxcofrb) = 0.0D+00
  grbsav(1:maxcofrb) = 0.0D+00
  grbsen(1:maxcofrb) = 0.0D+00
  grbtay(1:maxcofrb) = 0.0D+00
  gridx = 'uniform'
  gridy = 'uniform'
  hx = 0.0D+00
  hy = 0.0D+00
  ibs = 0
  ibump = 0
  icolrb(1:maxcofrb) = 0
  ierror = 0
  ifs = 0
  ihi = 0
  ijac = 1
  ilo = 0
  indx(1:3,1:maxnp) = 0
  iopt(1:maxpar) = 0
  ipar = 0
  ipivfl(1:maxnfl) = 0
  ipivrb(1:maxcofrb) = 0
  isotri(1:maxelm) = 0
  iwrite = 0
  jhi = 0
  jlo = 0
  maxnew = 10
  maxopt = 0
  maxsim = 10
  nbcrb = 0
  ncofrb = 0
  nelem = 0
  neqnfl = 0
  nferb = 0
  nlband = 0
  node(1:6,1:maxelm) = 0
  nodelm(1:maxnp) = 0
  np = 0
  npar = 1
  nparb = 0
  nparf = 0
  npe = 0

  do i = 1, 2*maxny-1
    nprof(i) = 0
  end do

  nsenfl = 5
  ntay = 0
  numnew = 0
  numopt = 0
  numsim = 0
  nx = 0
  ny = 0
  par(1:maxpar) = 0.0D+00
  parafl(1:maxpar) = 0.0D+00
  pararb(1:maxpar) = 0.0D+00
  pardif (1:maxpar) = 0.0D+00
  parrb(1:maxpar) = 0.0D+00
  parsav(1:maxpar) = 0.0D+00
  parsen(1:maxpar) = 0.0D+00
  partar(1:maxpar) = 0.0D+00

  do i = 1, 3
    do j = 1, 6
      do k = 1, 10
        phifl(i,j,k,1:maxelm) = 0.0D+00
      end do
    end do
  end do

  do i = 1,3
    do j = 1,maxcofrb
      do k = 1,15
        phirb(i,j,k,1:maxelm) = 0.0D+00
      end do
    end do
  end do

  do i = 1,maxnfl
    do j = 1,maxcofrb
      if ( i == j) then
        rb(i,j) = 1.0D+00
      else
        rb(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1,maxcofrb
    do j = 1,maxcofrb
      if ( i == j) then
        rbase(i,j) = 1.0D+00
      else
        rbase(i,j) = 0.0D+00
      end if
    end do
  end do

  region = ' '
  resfl(1:maxnfl) = 0.0D+00
  resflsav(1:maxnfl) = 0.0D+00
  resrb(1:maxcofrb) = 0.0D+00
  reynld = 1.0D+00
  reytay = 1.0D+00

  do i = 1,maxnfl
    do j = 1,maxcofrb
      if ( i == j) then
        senfl(i,j) = 1.0D+00
      else
        senfl(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1,maxcofrb
    do j = 1,maxcofrb
      if ( i == j) then
        senrb(i,j) = 1.0D+00
      else
        senrb(i,j) = 0.0D+00
      end if
    end do
  end do

  do i = 1,maxparb+2
    splbmp(i) = 0.0D+00
  end do

  splflo(1:maxparf) = 0.0D+00
  taubmp(1:maxparb+2) = 0.0D+00
  tauflo(1:maxparf) = 0.0D+00
  tecfil = 'tecplot.dat'
  tolnew = 0.0D+00
  tolopt = 0.0D+00
  tolsim = 0.0D+00
  value = 0.0D+00
  wateb = 0.0D+00
  watep = 0.0D+00
  wateu = 0.0D+00
  watev = 0.0D+00
  wquad(1:3) = 0.0D+00
  xbl = 0.0D+00
  xbr = 0.0D+00
  xc(1:maxnp) = 0.0D+00
  xmax = 0.0D+00
  xmin = 0.0D+00
  xprof = 0.0D+00
  xquad(1:3,1:maxelm) = 0.0D+00
  xrange = 0.0D+00
  xsiq(1:3) = 0.0D+00
  ybl = 0.0D+00
  ybr = 0.0D+00
  yc(1:maxnp) = 0.0D+00
  ymax = 0.0D+00
  ymin = 0.0D+00
  yquad(1:3,1:maxelm) = 0.0D+00
  yrange = 0.0D+00

  return
end
subroutine newtfl ( afl, area, eqn, gfl, gflafl, ierror, ifs, ijac,&
  indx, ipivfl, iwrite, ldafl, maxelm, maxnew, nelem, neqnfl, nlband, &
  node, np, npar, nparf, numnew, par, parafl, phifl, &
  region, resfl, rmax, splflo, tauflo, tolnew, xrange, yc, yrange )

!*****************************************************************************80
!
!! NEWTFL applies Newton's method to solve the full system.
!
!  Discussion:
!
!    The exact solution would have a zero residual, as computed by
!    the routine FXFL.  NEWTFL uses Newton's method to seek a solution
!    whose maximum residual is no more than TOLNEW.  The routine FPFL
!    is used to compute the Jacobian of the residual functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 8 ) AFL(LDAFL,NEQNFL),
!    AFL contains the Jacobian matrix for the full system,
!    stored in LINPACK general band storage mode.
!    The two dimensional array is of logical dimensions LDAFL by
!    NEQNFL.
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!  EQN    Input, character ( len = 2 ) EQN(NEQNFL).
!         EQN records the "type" of each equation that will be generated, and
!         which is associated with an unknown.  Note that most boundary
!         conditions do not result in an equation.  The current values are:
!
!         'U'  The horizontal momentum equation.
!         'UB' The condition U = 0 applied at a node on the bump.
!         'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!         'UW' The condition U = 0 applied at a node on a fixed wall.
!
!         'V'  The vertical momentum equation.
!         'VB' The condition V = 0 applied at a node on the bump.
!         'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!         'VW' The condition V = 0 applied at a node on a fixed wall.
!
!         'P'  The continuity equation.
!         'PB' The condition P = 0 applied at (XMAX,YMAX).
!
!  GFL    Input/output, real ( kind = 8 ) GFL(NEQNFL), the current solution
!         estimate for the full problem.
!
!  IERROR Output, integer ( kind = 4 ) IERROR, error flag.
!         0, no error occurred.
!         1, an error occurred, and the improved solution could not
!         be computed.
!
!  IFS    Input, integer ( kind = 4 ) IFS.
!         1, the inflow is modeled by C0 linear splines.
!         2, the inflow is modeled by C0 quadratic splines.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  IPIVFL Workspace, integer IPIVFL(NEQNFL), pivot vector for the solution
!         of the full linear system.
!
!  LDAFL  Input, integer ( kind = 4 ) LDAFL, the first dimension of the matrix AFL.
!
!  MAXNEW Input, integer ( kind = 4 ) MAXNEW, the maximum number of Newton steps
!         that may be taken.  10 should usually be enough.
!
!  MAXNFL Input, integer ( kind = 4 ) MAXNFL.
!
!         The maximum number of equations allowed for the full system.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of equations in the full system.
!
!  NLBAND Input, integer ( kind = 4 ) NLBAND.
!
!         The lower bandwidth of the matrix AFL.  The zero structure of AFL
!         is assumed to be symmetric, and so NLBAND is also the upper
!         bandwidth of AFL. 
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM).
!
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  NPAR   Input, integer ( kind = 4 ) NPAR.
!
!         The number of parameters.  NPAR = NPARF + NPARB + 1.
!
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!
!         The number of parameters associated with the position and
!         shape of the bump.
!
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1. 
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!
!         PAR is the current estimate for the parameters.
!
!  PHIFL  Input, real ( kind = 8 ) PHIFL(3,6,10,NELEM). 
!
!         PHIFL contains the value of a finite element basis function, its
!         derivative, or other information, evaluated at the quadrature
!         points.
!
!         The meaning of the entry PHIFL(I,J,K,L) is as follows. 
!         For the quadrature point I, and basis function J, in element L,
!         PHIFL(I,J,K,L) represents the value of:
!
!           K =  1, W, the finite element basis function for velocities;
!           K =  2, dWdX, the X derivative of W;
!           K =  3, dWdY, the Y derivative of W;
!           K =  4, Q, the finite element basis function for pressures;
!           K =  5, dQdX, the X derivative of Q;
!           K =  6, dQdY, the Y derivative of Q;
!           K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!           K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!           K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!           K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!         In particular, PHIFL(I,J,K,L) is the value of the quadratic
!         basis function W associated with local node J in element L,
!         evaluated at quadrature point I.
!
!         Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!         since there are only three linear basis functions.
!
!  REGION
!    Input, character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottome, with tangential velocity specifications
!    there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!  RESFL  Workspace, real ( kind = 8 ) RESFL(NEQNFL), the residual in the
!         full basis equations.
!
!  SPLFLO Input, real ( kind = 8 ) SPLFLO(NPARF).
!         SPLFLO contains the spline coefficients for the inflow.
!
!    Input, real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow.  There are NPARF of them, because the end
!    values of the spline are constrained to have particular
!    values.
!
!    Input, real ( kind = 8 ) TOLNEW, the Newton tolerance.
!    NEWTFL is asked to find an approximate solution so that
!    the maximum absolute value of all the residuals is no more
!    than TOLNEW.  A value such as 10E-7 is often reasonable,
!    though this depends on the actual equations being solved.
!
!    Input, real ( kind = 8 ) XC(NP).
!    The X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YC(NP).
!    The Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) afl(ldafl,neqnfl)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dmax
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gflafl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) idmax
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivfl(neqnfl)
  integer ( kind = 4 ) irmax
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) ixmax
  logical lmat
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) numnew
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) parafl(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmax0
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax0
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  Force the jacobian matrix to be evaluated on the first iteration.
!
  lmat = .false.
!
!  If the first Newton iteration failed, you may want to try again
!  by coming back here.
!
10    continue

  ierror = 0
  numnew = 0
!
!  Compute the norm of the initial solution estimate.
!
  ixmax = idamax(neqnfl,gfl,1)
  xmax = abs(gfl(ixmax))
  xmax0 = xmax
!
!  Evaluate the residual of the initial solution.
!
  call fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar, &
    nparf,par,phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)

  irmax = idamax(neqnfl,resfl,1)
  rmax = abs(resfl(irmax))
  rmax0 = rmax

  if ( 2 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' Step   MxNorm(X)   IXmax  MxNorm(FX)   IRmax'
    write(*,'(i6,g14.6,i6,g14.6,i6)')numnew,xmax,ixmax,rmax,irmax
  end if
!
!  Begin the Newton iteration.
!
  do numnew = 1,maxnew
!
!  If we have a valid, factored jacobian already, then we may
!  reuse it, if it's not too old, and if we're allowed.
!
    if ( 1 < ijac ) then
      if ( mod(numnew-1,ijac) == 0) then
        lmat = .false.
      else
        lmat = .true.
      end if
    else
      lmat = .false.
    end if
!
!  If it's time, evaluate and factor the jacobian.
!
    if ( .not. lmat ) then

      parafl(1:npar) = par(1:npar)
      gflafl(1:neqnfl) = gfl(1:neqnfl)

      call fpfl(afl,area,eqn,gfl,indx,ldafl,maxelm,nelem,neqnfl, &
        nlband,node,np,npar,par,phifl)

      call dfacfl(afl,ldafl,neqnfl,nlband,nlband,ipivfl,info)

      if ( info /= 0) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'NEWTFL - Fatal error!'
        write ( *, '(a)' ) '  The jacobian is singular.'
        write ( *, '(a,i6)' ) '  DFACFL returns INFO = ',info
        ierror = 1
        return
      else
        lmat = .true.
      end if

    end if
!
!  Solve the linear system A*DX = RES
!
    call dsolfl(afl,ldafl,neqnfl,nlband,nlband,ipivfl,resfl)

    idmax = idamax(neqnfl,resfl,1)
    dmax = abs(resfl(idmax))
!
!  Update the estimated solution G.
!
    gfl(1:neqnfl) = gfl(1:neqnfl) - resfl(1:neqnfl)
!
!  Compute the norm of the current solution.
!
    ixmax = idamax(neqnfl,gfl,1)
    xmax = abs(gfl(ixmax))
!
!  Evaluate the residual of the current estimated solution.
!
    call fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar, &
      nparf,par,phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)

    irmax = idamax(neqnfl,resfl,1)
    rmax = abs(resfl(irmax))

    if ( 2 <= iwrite ) then
      write(*,'(i6,g14.6,i6,g14.6,i6)')numnew,xmax,ixmax,rmax,irmax
    end if
!
!  Accept the iterate if the residual is small enough.
!
    if ( rmax <= tolnew) then
      return
    end if
!
!  Reject the iterate if the residual has grown too large.
!
    if ( 10.0D+00 *(rmax0+tolnew) < rmax .and. 1 < numnew ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEWTFL - Warning!'
      write ( *, * ) '  Residual too big on step ',numnew
      write ( *, * ) '  MxNorm of first FX = ',rmax0
      write ( *, * ) '  MxNorm of this FX =  ',rmax
      go to 20
    end if

  end do
!
!  The iteration has failed to converge, or may actually
!  have been terminated early.
!
20    continue

  ierror = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEWTFL - Warning!'
  write ( *, * ) '  No Newton convergence after   ',maxnew,' steps.'
  write ( *, * ) '  MxNorm of last step = ',dmax
  write ( *, * ) '  MxNorm of first X = ',xmax0
  write ( *, * ) '  MxNorm of last X =  ',xmax
  write ( *, * ) '  MxNorm of first FX  = ',rmax0
  write ( *, * ) '  MxNorm of last FX = ',rmax
  write ( *, * ) '  Tolerance for FX =    ',tolnew

  return
end
subroutine newtrb(arb,area,grb,grbarb,ierror,indx,ipivrb, &
  iwrite,maxcofrb,maxelm,maxnew,maxnfl,nbcrb,ncofrb,nelem, &
  nferb,node,np,npar,nparf,nx,ny,par,pararb,phirb, &
  rb,resrb,rmax,tauflo,tolnew,xc,xrange,yc,yrange)

!*****************************************************************************80
!
!! NEWTRB applies the Newton method to the reduced nonlinear state equations.
!
!  Discussion:
!
!    The routine is given an initial estimate of the solution of the reduced
!    nonlinear state equations, and seeks a better solution.
!
!    The exact solution would have a zero residual, as computed by
!    the routine FXRB.  NEWTRB uses Newton's method to seek a solution
!    whose maximum residual is no more than TOLNEW.  The routine FPRB
!    is used to compute the Jacobian of the residual functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 August 1996.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 8 ) ARB(MAXNRB,MAXNRB).
!    ARB contains the Jacobian matrix for the reduced basis system.
!
!    Input/output, real ( kind = 8 ) GRB(NCOFRB), the current solution
!    estimate for the reduced basis problem.
!
!    Output, integer ( kind = 4 ) IERROR, error flag.
!    0, no error occurred.
!    1, an error occurred, and the improved solution could not be computed.
!
!    Workspace, integer IPIVRB(NCOFRB), pivot vector for the solution
!    of the reduced linear system.
!
!    Input, integer ( kind = 4 ) IWRITE.
!    IWRITE controls the amount of output printed.
!    0 = little, 1=some, 2=a lot.
!
!    Input, integer ( kind = 4 ) MAXNEW, the maximum number of Newton steps
!    that may be taken.  10 should usually be enough.
!
!    Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!    Input, integer ( kind = 4 ) NCOFRB, the number of equations in the reduced system.
!
!    Input, integer ( kind = 4 ) NPAR.
!    The number of parameters.  NPAR = NPARF + NPARB + 1.
!    The parameters control the shape of the inflow,
!    the shape of the bump obstacle, and the strength of the flow.
!
!    Input, real ( kind = 8 ) PAR(NPAR).
!    PAR is the current estimate for the parameters.
!
!    Output, real ( kind = 8 ) PARMAT(NPAR).
!    PARMAT contains the parameters where the Jacobian was generated.
!
!    Input, real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!    PHIRB contains the values of a finite element basis function
!    or its X or Y derivative, in a given element, at a given
!    quadrature point, for a particular reduced basis function.
!
!    For PHIRB(I,J,K,L), index J refers to the reduced basis
!    basis functions, for J = 0 to NCOFRB.
!
!    The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!    For the quadrature point I, and reduced basis function J,
!    in element L, PHIRB(I,J,K,L) represents the value of:
!
!    K = 1, WUrb, the finite element U velocity basis function;
!    K = 2, dWUrbdX, the X derivative of WUrb;
!    K = 3, dWUrbdY, the Y derivative of WUrb;
!    K = 4, WVrb, the finite element V velocity basis function;
!    K = 5, dWVrbdX, the X derivative of WVrb;
!    K = 6, dWVrbdY, the Y derivative of WVrb;
!    K = 7, Q, the finite element pressure basis function.
!    K = 8, dQrbdX, the X derivative of Qrb;
!    K = 9, dQrbdY, the Y derivative of Qrb.
!    K = 10, WU0rb, same as WUrb, with zero BC.
!    K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!    K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!    K = 13, WV0rb, same as WVrb, with zero BC.
!    K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!    K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    Workspace, real ( kind = 8 ) RESRB(NCOFRB), the residual in the
!    reduced basis equations, evaluated at the coefficient
!    vector GRB.
!
!    Input, real ( kind = 8 ) TOLNEW, the Newton tolerance.
!    NEWTRB is asked to find an approximate solution so that
!    the maximum absolute value of all the residuals is no more
!    than TOLNEW.  A value such as 10E-7 is often reasonable,
!    though this depends on the actual equations being solved.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dmax
  real ( kind = 8 ) grb(ncofrb)
  real ( kind = 8 ) grbarb(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) idmax
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivrb(ncofrb)
  integer ( kind = 4 ) irmax
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) numnew
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) pararb(npar)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmax0
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax0
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
!  If the first Newton iteration failed, you may want to try again
!  by coming back here.
!
10    continue

  ierror = 0
  numnew = 0
!
!  Compute the norm of the initial solution estimate.
!
  ixmax = idamax(ncofrb,grb,1)
  xmax = abs(grb(ixmax))
  xmax0 = xmax
!
!  Evaluate the residual of the initial solution.
!
  reynld = par(npar)

  call fxrb(area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
    nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb, &
    resrb,reynld,tauflo,xc,xrange,yc,yrange)

  irmax = idamax(ncofrb,resrb,1)
  rmax = abs(resrb(irmax))
  rmax0 = rmax

  if ( 2 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' Step   MxNorm(X)   IXmax  MxNorm(FX)   IRmax'
    write(*,'(i6,g14.6,i6,g14.6,i6)')numnew,xmax,ixmax,rmax,irmax
  end if
!
!  Begin the Newton iteration.
!
  do numnew = 1,maxnew
!
!  Evaluate the Jacobian.
!
    pararb(1:npar) = par(1:npar)
    grbarb(1:ncofrb) = grb(1:ncofrb)

    call fprb(arb,area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb, &
      ncofrb,nelem,nferb,node,np,nx,ny,phirb,rb,reynld,xc,xrange,yc,yrange)
!
!  Factor the Jacobian.
!
    call dfacrb(arb,maxcofrb,ncofrb,ipivrb,info)

    if ( info /= 0) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEWTRB - Fatal error!'
      write ( *, '(a)' ) '  The reduced Jacobian is singular.'
      write ( *, * ) '  DFACRB returns INFO = ',info
      ierror = 1
      return
    end if
!
!  Solve the linear system A*DX = RES
!
    call dsolrb(arb,maxcofrb,ncofrb,ipivrb,resrb)

    idmax = idamax(ncofrb,resrb,1)
    dmax = abs(resrb(idmax))
!
!  Update the estimated solution G.
!
    grb(1:ncofrb) = grb(1:ncofrb) - resrb(1:ncofrb)
!
!  Compute the norm of the current solution.
!
    ixmax = idamax(ncofrb,grb,1)
    xmax = abs(grb(ixmax))
!
!  Evaluate the residual of the current estimated solution.
!
    call fxrb(area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
      nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb, &
      resrb,reynld,tauflo,xc,xrange,yc,yrange)

    irmax = idamax(ncofrb,resrb,1)
    rmax = abs(resrb(irmax))

    if ( 2 <= iwrite ) then
      write(*,'(i6,g14.6,i6,g14.6,i6)')numnew,xmax,ixmax,rmax,irmax
    end if
!
!  Accept the iterate if the residual is small enough.
!
    if ( rmax <= tolnew) then
      return
    end if
!
!  Reject the iterate if the residual has grown too large.
!
    if ( 10.0D+00 *(rmax0+tolnew) < rmax .and. 1 < numnew ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEWTRB - Warning!'
      write ( *, * ) '  Residual too big on step ',numnew
      write ( *, * ) '  MxNorm of first FX = ',rmax0
      write ( *, * ) '  MxNorm of this FX =  ',rmax
      go to 20
    end if

  end do
!
!  The iteration has failed to converge, or may actually
!  have been terminated early.
!
20    continue

  ierror = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEWTRB - Warning!'
  write ( *, * ) '  No Newton convergence after   ',numnew,' steps.'
  write ( *, '(a)' ) ' '
  write ( *, * ) '  MxNorm of first X = ',xmax0
  write ( *, * ) '  MxNorm of first FX  = ',rmax0
  write ( *, '(a)' ) ' '
  write ( *, * ) '  MxNorm of last step = ',dmax
  write ( *, * ) '  MxNorm of last X =  ',xmax
  write ( *, * ) '  MxNorm of last FX = ',rmax

  return
end
subroutine optdiffl(afl,area,cost,dopt,eqn,etaq,gfl,gflafl, &
  gflopt,gfltar,gridx,gridy,ibs,ierror,ifs,ijac,indx,iopt, &
  ipivfl,isotri,ivopt,iwrite,ldafl,liv,lv,maxelm,maxnew,maxnfl, &
  maxnp,maxny,maxopt,maxpar,maxparb,maxparf,maxsim,nelem, &
  neqnfl,nlband,node,nodelm,np,npar,nparb,nparf,nprof,numdif, &
  numopt,nx,ny,par,parafl,paropt,phifl,region,resfl,splbmp, &
  splflo,taubmp,tauflo,tolnew,tolopt,tolsim,vopt,wateb,watep, &
  wateu,watev,wquad,xbl,xbr,xc,xopt,xquad,xrange,xsiq,ybl, &
  ybr,yc,yquad,yrange)

!*****************************************************************************80
!
!! OPTDIFFL optimizes the full problem, without gradient information.
!
!  Discussion:
!
!    OPTDIFFL searches for a set of parameters PAROPT,
!    and the corresponding flow solution GFLOPT, which minimize
!    the cost function COST.
!
!    The ACM TOMS 611 routine SNOIT is used, which does not require
!    direct information about the gradient of COST with respect to
!    the parameters PAROPT.  Instead, it estimates this information
!    indirectly, via finite differences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) AFL(LDAFL,MAXNFL).
!    If Newton iteration is being carried out, AFL contains the
!    Jacobian matrix for the full system.
!    If Picard iteration is being carried out, AFL contains the
!    Picard matrix for the full system.
!    AFL is stored in LINPACK general band storage mode, with
!    logical dimensions (3*NBANDL+1, NEQNFL).
!
!    Input, real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Output, real ( kind = 8 ) COST.
!    COST contains the current value of the cost function.  This
!    is the function which the optimizer is to minimize.
!    COST = WATEP*COSTP + WATEB*COSTB + WATEU*COSTU + WATEV*COSTV
!
!    Workspace, real ( kind = 8 ) DOPT(MAXPAR).
!    DOPT contains scaling factors used during an optimization.
!    These scaling factors are intended to adjust problems
!    in which some variables are typically very much smaller
!    or larger than others.
!
!  EQN    Input, character ( len = 2 ) EQN(MAXNFL).
!         EQN records the "type" of each equation that will be generated, and
!         which is associated with an unknown.  Note that most boundary
!         conditions do not result in an equation.  The current values are:
!
!         'U'  The horizontal momentum equation.
!         'UB' The condition U = 0 applied at a node on the bump.
!         'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!         'UW' The condition U = 0 applied at a node on a fixed wall.
!         'U0' A dummy value of U = 0 should be set.
!
!         'V'  The vertical momentum equation.
!         'VB' The condition V = 0 applied at a node on the bump.
!         'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!         'VW' The condition V = 0 applied at a node on a fixed wall.
!         'V0' A dummy value of V = 0 should be set.
!
!         'P'  The continuity equation.
!         'PB' The condition P = 0 applied at (XMAX,YMAX).
!         'P0' A dummy value of P = 0 should be set.
!
!  ETAQ   Input, real ( kind = 8 ) ETAQ(3).
!         ETAQ contains the "Eta" coordinates of the quadrature points.
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!         GFL contains the current solution estimate for the full problem,
!         containing the pressure and velocity coefficients.
!         The vector INDX must be used to index this data.
!
!  GFLAFL Output, real ( kind = 8 ) GFLAFL(NEQNFL).
!         GFLAFL stores the value of GFL at which the Jacobian
!         was generated.
!
!  GFLOPT Output, real ( kind = 8 ) GFLOPT(NEQNFL).
!         GFLOPT stores the value of a full solution which is being
!         optimized.
!
!  GFLTAR Input, real ( kind = 8 ) GFLTAR(NEQNFL).
!         GFLTAR is a target solution, used to generate data that defines
!         the cost functional.  The corresponding parameters are PARTAR.
!
!  IBS    Input, integer ( kind = 4 ) IBS.
!         IBS is the bump shape option.
!         0, piecewise constant function.
!         1, piecewise linear function.
!         2, piecewise quadratic function.
!
!  IERROR Output, integer ( kind = 4 ) IERROR.
!         0, the optimization was successful.
!         1, the optimization failed.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP).
!         INDX(I,J) contains, for each node J, the global index of U,
!         V and P at that node, or 0 or a negative value.  The global
!         index of U, V, or P is the index of the coefficient vector
!         that contains the value of the finite element coefficient
!         associated with the corresponding basis function at the
!         given node.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  IOPT   Workspace, integer IOPT(MAXPAR).
!         IOPT is used during an optimization.  For each parameter I,
!         the meaning of IOPT(I) is:
!         0, the parameter value must remain fixed;
!         1, the parameter value may be varied.
!
!  IPIVFL Workspace, integer IPIVFL(NEQNFL).
!         IPIVFL is a pivot vector for the solution of the full
!         linear system.
!
!  ISOTRI Input, integer ( kind = 4 ) ISOTRI(NELEM).
!         0, the element is NOT isoparametric, and the nodes never move.
!         That means that the quadrature points are only computed once.
!
!         1, the element is NOT isoparametric, but the nodes may move.
!         Quadrature point locations must be updated on each step.
!         This could occur for elements above, but not touching, the bump.
!
!         2, the element is isoparametric.
!
!  IVOPT  Workspace, integer IVOPT(LIV).
!         IVOPT provides integer workspace for several of the
!         optimization routines.
!
!  IWRITE Input, integer ( kind = 4 ) IWRITE.
!         IWRITE controls the amount of output printed.
!         0, print out the least amount.
!         1, print out some.
!         2, print out a lot.
!
!  LDAFL  Input, integer ( kind = 4 ) LDAFL.
!         LDAFL is the first dimension of the matrix AFL as declared in
!         the main program.  LDAFL must be at least 3*NLBAND+1.
!
!  LIV    Input, integer ( kind = 4 ) LIV.
!         LIV is the dimension of the work vector IVOPT, used by
!         the ACM TOMS 611 optimization package.  LIV is always 60.
!
!  LV     Input, integer ( kind = 4 ) LV.
!         LV is the dimension of the work vector VOPT, used by
!         the ACM TOMS 611 optimization package. 
!
!  MAXELM Input, integer ( kind = 4 ) MAXELM.
!         MAXELM is the maximum number of elements.
!
!  MAXNEW Input, integer ( kind = 4 ) MAXNEW.
!         MAXNEW is the maximum number of steps to take in one Newton
!         iteration.  A typical value is 20.
!
!  MAXNFL Input, integer ( kind = 4 ) MAXNFL.
!         MAXNFL is the maximum number of equations or coefficients allowed
!         for the full system.  MAXNFL must be used instead of NEQNFL as
!         the leading dimension of certain multi-dimensional arrays.
!
!  MAXNP  Input, integer ( kind = 4 ) MAXNP.
!         MAXNP is the maximum number of nodes allowed in the program.
!
!  MAXNY  Input, integer ( kind = 4 ) MAXNY.
!         MAXNY is the maximum size of NY that the program can handle.
!
!  MAXOPT Input, integer ( kind = 4 ) MAXOPT.
!         MAXOPT is the maximum number of optimization steps.
!
!  MAXPAR Input, integer ( kind = 4 ) MAXPAR.
!         MAXPAR is the maximum number of parameters allowed.
!         MAXPAR = MAXPARF + MAXPARB + 1.
!
!  MAXPARB
!         Input, integer ( kind = 4 ) MAXPARB.
!         MAXPARB is the maximum number of bump parameters allowed.
!
!  MAXPARF
!         Input, integer ( kind = 4 ) MAXPARF.
!         MAXPARF is the maximum number of inflow parameters allowed.
!
!  MAXSIM Input, integer ( kind = 4 ) MAXSIM.
!         MAXSIM is the maximum number of steps to take in one Picard
!         iteration.  A typical value is 20.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM.
!         NELEM is the number of elements.
!         NELEM can be determined as 2*(NX-1)*(NY-1).
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL.
!         NEQNFL is the number of equations (and coefficients) in the full
!         finite element system.
!
!  NLBAND Input, integer ( kind = 4 ) NLBAND.
!         NLBAND is the lower bandwidth of the matrix AFL.
!         The zero structure of AFL is assumed to be symmetric, and so
!         NLBAND is also the upper bandwidth of AFL.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!         NODE(I,J) contains, for an element J, the global index of
!         the node whose local number in J is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!           Global nodes   Elements      NODE
!                                                          1  2  3  4  5  6
!           74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                          |  /   /|
!           73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                          |/   /  |
!           72  82  92     2   1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP.
!         NP is the number of nodes used to define the finite element mesh.
!         Typically, the mesh is generated as a rectangular array, with
!         an odd number of nodes in the horizontal and vertical directions.
!         The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!  NPAR   Input, integer ( kind = 4 ) NPAR.
!         NPAR is the number of parameters.
!
!         NPAR = NPARF + NPARB + 1.
!
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!         NPARB is the number of parameters associated with the position and
!         shape of the bump.
!
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1.
!
!  NPROF  Input, integer ( kind = 4 ) NPROF(2*MAXNY-1).
!         NPROF contains the numbers of the nodes along the profile
!         line.
!
!  NX     Input, integer ( kind = 4 ) NX.
!         NX controls the spacing of nodes and elements in
!         the X direction.  There are 2*NX-1 nodes along various
!         lines in the X direction.
!
!         Roughly speaking, NX (or 2*NX) is the number of elements along
!         a line in the X direction.
!
!  NY     Input, integer ( kind = 4 ) NY.
!         NY controls the spacing of nodes and elements in
!         the Y direction.  There are 2*NY-1 nodes along various
!         lines in the Y direction.
!
!         Roughly speaking, NY (or 2*NY) is the number of elements along
!         a line in the Y direction.
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!         PAR is the current estimate for the parameters.
!
!         PAR(1:NPARF)             = inflow controls.
!
!         PAR(NPARF+1:NPARF+NPARB) = bump controls.
!
!         PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!  PARAFL Output, real ( kind = 8 ) PARAFL(NPAR).
!         PARAFL contains the parameters where the Picard matrix or
!         Jacobian of the full system was generated.
!
!  PAROPT Output, real ( kind = 8 ) PAROPT(NPAR).
!         PAROPT contains the estimate for the optimizing parameter
!         values which minimize the cost.
!
!  PHIFL  Input, real ( kind = 8 ) PHIFL(3,6,10,NELEM).
!         PHIFL contains the value of a finite element basis function, its
!         derivative, or other information, evaluated at the quadrature
!         points.
!
!         The meaning of the entry PHIFL(I,J,K,L) is as follows.
!         For the quadrature point I, and basis function J, in element L,
!         PHIFL(I,J,K,L) represents the value of:
!
!           K =  1, W, the finite element basis function for velocities;
!           K =  2, dWdX, the X derivative of W;
!           K =  3, dWdY, the Y derivative of W;
!           K =  4, Q, the finite element basis function for pressures;
!           K =  5, dQdX, the X derivative of Q;
!           K =  6, dQdY, the Y derivative of Q;
!           K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!           K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!           K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!           K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!         In particular, PHIFL(I,J,K,L) is the value of the quadratic
!         basis function W associated with local node J in element L,
!         evaluated at quadrature point I.
!
!         Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!         since there are only three linear basis functions.
!
!    Input, character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottome, with tangential velocity specifications
!    there.
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    Workspace, real ( kind = 8 ) RESFL(NEQNFL).
!    RESFL contains the residual in the full basis equations.
!
!    Input, real ( kind = 8 ) SPLBMP(NPARB+2).
!    SPLBMP contains the spline coefficients for the bump.
!
!    Input, real ( kind = 8 ) SPLFLO(NPARF).
!    SPLFLO contains the spline coefficients for the inflow.
!
!  TAUBMP Input, real ( kind = 8 ) TAUBMP(NPARB+2).
!         TAUBMP contains the location of the spline abscissas for
!         the bump.  There are NPARB+2 of them, because the end values
!         of the spline are constrained to have particular values.
!
!  TAUFLO Input, real ( kind = 8 ) TAUFLO(NPARF).
!         TAUFLO contains the location of the spline abscissas for
!         the inflow.  There are NPARF of them, because the end
!         values of the spline are constrained to have particular
!         values.
!
!  TOLNEW Input, real ( kind = 8 ) TOLNEW.
!         TOLNEW is the convergence tolerance for the Newton
!         iteration.
!
!  TOLOPT Input, real ( kind = 8 ) TOLOPT.
!         TOLOPT is the convergence tolerance for the optimization.
!         If TOLOPT is zero, then default values are used.
!
!    Input, real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the Picard
!    iteration.
!
!    Workspace, real ( kind = 8 ) VOPT(LV).
!    VOPT provides real workspace for the optimization routines.
!
!    Input, real ( kind = 8 ) WATEB.
!    WATEB is the multiplier of the bump control cost used
!    when computing the total cost.
!
!    Input, real ( kind = 8 ) WATEP, WATEU, WATEV.
!    WATEP, WATEU and WATEV are weights used in computing the
!    cost function based on the costs of the flow discrepancy.
!
!    Input, real ( kind = 8 ) WQUAD(3).
!    WQUAD contains the weights for Gaussian quadrature.
!
!    Input, real ( kind = 8 ) XBL.
!    XBL is the X coordinate of the left corner of the bump.
!
!    Input, real ( kind = 8 ) XBR.
!    XBR is the X coordinate of the right corner of the bump.
!
!    Input, real ( kind = 8 ) XC(NP).
!    XC contains the X coordinates of the nodes.
!
!    Workspace, real ( kind = 8 ) XOPT(MAXPAR).
!    XOPT is used by the optimization routines to hold only
!    the values of parameters which are allowed to vary.
!
!    Input, real ( kind = 8 ) XQUAD(3,NELEM).
!    The X coordinates of the quadrature points for
!    each element.
!
!    Input, real ( kind = 8 ) XRANGE.
!    The total width of the region.
!
!    Input, real ( kind = 8 ) XSIQ(3).
!    The "Xsi" coordinates of the quadrature points.
!
!    Input, real ( kind = 8 ) YBL.
!    The Y coordinate of the left corner of the bump.
!
!    Input, real ( kind = 8 ) YBR.
!    YBR is the Y coordinate of the right corner of the bump.
!
!    Input, real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YQUAD(3,NELEM).
!    The Y coordinates of the quadrature points for
!    each element.
!
!    Input, real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
!  Set parameters that are independent.
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) liv
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) maxelm   
  integer ( kind = 4 ) maxnfl  
  integer ( kind = 4 ) maxnp              
  integer ( kind = 4 ) maxny
  integer ( kind = 4 ) maxpar     
  integer ( kind = 4 ) maxparb
  integer ( kind = 4 ) maxparf
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) afl(ldafl,maxnfl)
  real ( kind = 8 ) area(3,maxelm)
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) dopt(maxpar)
  character ( len = 2 ) eqn(maxnfl)
  real ( kind = 8 ) etaq(3)
  real ( kind = 8 ) gfl(maxnfl)
  real ( kind = 8 ) gflafl(maxnfl)
  real ( kind = 8 ) gflopt(maxnfl)
  real ( kind = 8 ) gfltar(maxnfl)
  character ( len = 20 ) gridx
  character ( len = 20 ) gridy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) indx(3,maxnp)
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) ipivfl(maxnfl)
  integer ( kind = 4 ) isotri(maxelm)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivopt(liv)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) maxsim
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,maxelm)
  integer ( kind = 4 ) nodelm(np)
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nprof(2*maxny-1)
  integer ( kind = 4 ) numdif
  integer ( kind = 4 ) numnew
  integer ( kind = 4 ) numopt
  integer ( kind = 4 ) numsim
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(maxpar)
  real ( kind = 8 ) parafl(maxpar)
  real ( kind = 8 ) paropt(maxpar)
  real ( kind = 8 ) phifl(3,6,10,maxelm)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(maxnfl)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) splbmp(maxparb+2)
  real ( kind = 8 ) splflo(maxparf)
  real ( kind = 8 ) taubmp(maxparb+2)
  real ( kind = 8 ) tauflo(maxparf)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) vopt(lv)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(maxnp)
  real ( kind = 8 ) xopt(maxpar)
  real ( kind = 8 ) xquad(3,maxelm)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yquad(3,maxelm)
  real ( kind = 8 ) yrange
!
  ierror = 0
!
!  Copy the initial solution estimate.
!
  paropt(1:npar) = par(1:npar)
  gflopt(1:neqnfl) = gfl(1:neqnfl)
!
!  Initialize the local optimization data.
!
  cost = 0.0D+00
  dopt(1:npar) = 1.0D+00
  ivopt(1:liv) = 0 
  nopt = 0
  vopt(1:lv) = 0.0D+00
  xopt(1:maxpar) = 0.0D+00
!
!  Set the TOMS 611 data to default values,
!  and then modify some values.
!
  ival = 2
  call deflt(ival,ivopt,liv,lv,vopt)

  vopt(31) = tolopt
  vopt(32) = tolopt
  vopt(33) = tolopt
  vopt(34) = tolopt
  vopt(37) = tolopt

  ivopt(1) = 12
  ivopt(19) = 0
!
!  Set the step counters.
!
  numdif = 0
  numopt = 0
!
!  Take the next optimization step.
!
  do

    if ( maxopt < numopt ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'OPTDIFFL - Warning!'
      write ( *, '(a)' ) '  The number of optimization steps was exceeded.'
      return
    end if
!
!  Make a "copy" of PAR that only contains the free variables.
!
    nopt = 0
    do i = 1,npar
      if ( iopt(i) == 1) then
        nopt = nopt+1
        xopt(nopt) = paropt(i)
      end if
    end do
!
!  Call the optimizer to get a new parameter estimate.
!
    call snoit(dopt,cost,ivopt,liv,lv,nopt,vopt,xopt)
!
!  Copy the new free variable values back into PAR.
!
    nopt = 0
    do i = 1,npar
      if ( iopt(i) == 1) then
        nopt = nopt+1
        paropt(i) = xopt(nopt)
      end if
    end do
!
!  For the given values of PAROPT, set up the flow problem.
!  We are only varying the REYNLD parameter, and no geometric
!  quantities vary with REYNLD, so we only have to make this
!  call once.
!
    if ( numopt == 0) then
      call setgeo(area,etaq,gridx,gridy,ibs,isotri,nelem,node, &
        nodelm,np,npar,nparb,nparf,nx,ny,paropt,phifl,region, &
        splbmp,taubmp,wquad,xbl,xbr,xc,xquad,xrange, &
        xsiq,ybl,ybr,yc,yquad,yrange)
    end if
!
!  Apply Picard's method to the approximate solution GFLOPT.
!
    call picfl(afl,area,eqn,gflopt,ierror,ifs,indx,ipivfl,iwrite,&
      ldafl,maxsim,nelem,neqnfl,nlband,node,np,npar,nparf, &
      numsim,paropt,phifl,region,resfl,rmax,splflo,tauflo, &
      tolsim,xc,xrange,yc,yrange)
!
!  Apply Newton's method to the approximate solution GFLOPT.
!
      if ( rmax <= tolnew) then
        write ( *, '(a)' ) 'OPTDIFFL - Picard iterate skips Newton.'
      else
        call newtfl(afl,area,eqn,gflopt,gflafl,ierror,ifs,ijac, &
          indx,ipivfl,iwrite,ldafl,maxelm,maxnew,nelem,neqnfl,nlband, &
          node,np,npar,nparf,numnew,paropt,parafl,phifl, &
          region,resfl,rmax,splflo,tauflo,tolnew,xrange,yc,yrange)

      if ( ierror /= 0) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'OPTDIFFL - Fatal error!'
        write ( *, '(a)' ) '  NEWTFL failed!'
        write ( *, '(a)' ) '  The parameters at which failure occurred:'
        write ( *, '(a)' ) ' '
        call prpar(iopt,npar,nparb,nparf,paropt)
        ierror = 1
        return
      end if

    end if
!
!  Compute the cost function COST.
!
    call getcst(cost,costb,costp,costu,costv,gflopt,gfltar, &
      indx,neqnfl,np,nparb,nprof,ny,splbmp, &
      taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)

    if ( ivopt(1) == 1 .and. 0 <= iwrite ) then
      call prpar(iopt,npar,nparb,nparf,paropt)
      write ( *, '(a,g14.6)' ) '  Cost = ',cost
    end if
!
!  If IVOPT(1) is 1, then this was a call for a legitimate
!  solution candidate.
!
!  If IVOPT(1) is 2, then this was a call for a temporary
!  solution used only for estimating the gradient.
!
!  Other values of IVOPT call for acceptance or rejection.
!
    if ( ivopt(1) == 1) then
      numopt = numopt+1
    else if ( ivopt(1) == 2) then
      numdif = numdif+1
    else if ( 3 <= ivopt(1) .and. ivopt(1) <= 8 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Convergence to a minimizer was achieved!'
      return
    else if ( 8 < ivopt(1) ) then
      write ( *, '(a)' ) ' '
      write ( *, * ) '  Bad value of IVOPT(1) = ',ivopt(1)
      return
    end if

  end do

end
subroutine optdifrb(arb,area,cost,dopt,gflrb,gfltar,gfltmp, &
  grb,grbarb,grbopt,ierror,indx,iopt,ipivrb, &
  ivopt,iwrite,liv,lv,maxcofrb,maxelm,maxnew,maxnfl,maxnp, &
  maxny,maxopt,maxpar,maxparb,maxsim,nbcrb,ncofrb,nelem,neqnfl, &
  nferb,node,np,npar,nparb,nparf,nprof,numdif,numopt,nx,ny,par, &
  pararb,paropt,phirb,rb,resrb,splbmp,tauflo,taubmp,tolnew, &
  tolopt,tolsim,vopt,wateb,watep,wateu,watev,xbl,xbr,xc,xopt, &
  xrange,ybl,ybr,yc,yrange)

!*****************************************************************************80
!
!! OPTDIFRB optimizes the reduced problem, without gradient information.
!
!  Discussion:
!
!    OPTDIFRB searches for a set of parameters PAROPT,
!    and the corresponding flow solution GRBOPT, which minimize
!    the cost function COST.
!
!    The ACM TOMS 611 routine SNOIT is used, which does not require
!    direct information about the gradient of COST with respect to
!    the parameters PAROPT.  Instead, it estimates this information
!    indirectly, via finite differences.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 8 ) ARB(MAXNRB,MAXNRB).
!    ARB contains the Jacobian or Picard matrix for the reduced
!    Navier Stokes system, stored as a dense NCOFRB by NCOFRB array.
!
!    Input, real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    real ( kind = 8 ) COST, the current value of the cost function.  This
!    is the function which the optimizer is to minimize.
!      COST = WATEP*COSTP + WATEB*COSTB + WATEU*COSTU + WATEV*COSTV
!
!    real ( kind = 8 ) DOPT(MAXPAR).
!    DOPT contains scaling factors used during an optimization.
!    These scaling factors are intended to adjust problems
!    in which some variables are typically very much smaller
!    or larger than others.
!
!    real ( kind = 8 ) GFLRB(NEQNFL).
!    GFLRB is the solution value at which the reduced basis was computed.
!    The corresponding parameters are PARRB.
!
!    real ( kind = 8 ) GFLTAR(NEQNFL).
!    GFLTAR is a target solution, used to generate data that defines
!    the cost functional.  The corresponding parameters are PARTAR.
!
!    Workspace, real ( kind = 8 ) GFLTMP(NEQNFL).
!
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    real ( kind = 8 ) GRBARB(NCOFRB).
!    GRBARB contains the reduced basis coefficients at which
!    the matrix ARB was last evaluated.
!
!    real ( kind = 8 ) GRBOPT(NCOFRB).
!    GRBOPT stores the value of a reduced solution which is being
!    optimized.
!
!    Workspace, real ( kind = 8 ) GRBTMP(NCOFRB).
!
!    integer ( kind = 4 ) IERROR.
!    IERROR is an error flag.
!    0, no error occurred in this routine.
!    nonzero, an error occurred.
!
!  INDX   integer ( kind = 4 ) INDX(3,NP).
!         INDX(I,J) contains, for each node J, the global index of U,
!         V and P at that node, or 0 or a negative value.  The global
!         index of U, V, or P is the index of the coefficient vector
!         that contains the value of the finite element coefficient
!         associated with the corresponding basis function at the
!         given node.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  IOPT   integer ( kind = 4 ) IOPT(MAXPAR).
!         IOPT is used during an optimization.  For each parameter I,
!         the meaning of IOPT(I) is:
!         0, the parameter value must remain fixed;
!         1, the parameter value may be varied.
!
!  IPIVRB Workspace, integer IPIVRB(NCOFRB).
!         IPIVRB is a pivot vector for the solution of the reduced
!         linear system.
!
!  IVOPT  integer ( kind = 4 ) IVOPT(LIV).
!         IVOPT provides integer workspace for several of the
!         optimization routines.
!
!  IWRITE integer IWRITE.
!         IWRITE controls the amount of output printed.
!         0, print out the least amount.
!         1, print out some.
!         2, print out a lot.
!
!  LIV    integer ( kind = 4 ) LIV.
!         LIV is the dimension of the work vector IVOPT, used by
!         the ACM TOMS 611 optimization package.  LIV is always 60.
!
!  LV     integer ( kind = 4 ) LV.
!         LV is the dimension of the work vector VOPT, used by
!         the ACM TOMS 611 optimization package. 
!
!  MAXELM integer MAXELM.
!         MAXELM is the maximum number of elements.
!
!  MAXNEW integer MAXNEW.
!         MAXNEW is the maximum number of steps to take in one Newton
!         iteration.  A typical value is 20.
!
!  MAXNFL integer MAXNFL.
!         MAXNFL is the maximum number of equations or coefficients allowed
!         for the full system.  MAXNFL must be used instead of NEQNFL as
!         the leading dimension of certain multi-dimensional arrays.
!
!  MAXNP  integer ( kind = 4 ) MAXNP.
!         MAXNP is the maximum number of nodes allowed in the program.
!
!  MAXNRB integer MAXNRB.
!         The maximum number of equations allowed for the reduced basis system.
!
!  MAXNY  integer ( kind = 4 ) MAXNY.
!         MAXNY is the maximum size of NY that the program can handle.
!
!  MAXOPT integer MAXOPT.
!         MAXOPT is the maximum number of optimization steps.
!
!  MAXPAR integer MAXPAR.
!         MAXPAR is the maximum number of parameters allowed.
!         MAXPAR = MAXPARF + MAXPARB + 1.
!
!  MAXPARB
!         integer ( kind = 4 ) MAXPARB.
!         MAXPARB is the maximum number of bump parameters allowed.
!
!  MAXSIM integer MAXSIM.
!         MAXSIM is the maximum number of steps to take in one Picard
!         iteration.  A typical value is 20.
!
!  NELEM  integer ( kind = 4 ) NELEM.
!         NELEM is the number of elements.
!         NELEM can be determined as 2*(NX-1)*(NY-1).
!
!  NEQNFL integer NEQNFL.
!         NEQNFL is the number of equations (and coefficients) in the full
!         finite element system.
!
!  NCOFRB integer NCOFRB.
!         NCOFRB is the number of basis functions, reduced state equations and
!         coefficients in the reduced basis system.
!
!  NP     integer ( kind = 4 ) NP.
!         NP is the number of nodes used to define the finite element mesh.
!         Typically, the mesh is generated as a rectangular array, with
!         an odd number of nodes in the horizontal and vertical directions.
!         The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!  NPAR   integer ( kind = 4 ) NPAR.
!         NPAR is the number of parameters.
!
!         NPAR = NPARF + NPARB + 1.
!
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  integer ( kind = 4 ) NPARB.
!         NPARB is the number of parameters associated with the position and
!         shape of the bump.
!
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  integer ( kind = 4 ) NPARF.
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1.
!
!  NPROF  integer ( kind = 4 ) NPROF(2*MAXNY-1).
!         NPROF contains the numbers of the nodes along the profile
!         line.
!
!  NUMDIF integer NUMDIF.
!         NUMDIF is the number of flow solutions generated strictly for
!         finite difference calculations.
!
!  NUMOPT integer NUMOPT.
!         NUMOPT is the number of flow solutions calculated during
!         an optimization which were actual candidate minimizers.
!
!  NY     integer ( kind = 4 ) NY.
!         NY controls the spacing of nodes and elements in
!         the Y direction.  There are 2*NY-1 nodes along various
!         lines in the Y direction.
!
!         Roughly speaking, NY (or 2*NY) is the number of elements along
!         a line in the Y direction.
!
!  PAR    real ( kind = 8 ) PAR(NPAR).
!         PAR is the current estimate for the parameters.
!
!         PAR(1:NPARF)             = inflow controls.
!
!         PAR(NPARF+1:NPARF+NPARB) = bump controls.
!
!         PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!  PARARB real ( kind = 8 ) PARARB(NPAR).
!         PARARB contains the parameters where the Picard matrix or
!         Jacobian of the reduced system was generated.
!
!  PAROPT real ( kind = 8 ) PAROPT(NPAR).
!         PAROPT contains the estimate for the optimizing parameter
!         values which minimize the cost.
!
!  PHIRB  Input, real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!  RB     real ( kind = 8 ) RB(MAXNFL,NCOFRB).
!         RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!
!         RB is generated by computing a finite element solution GFL.
!         A copy of this solution will be saved and called "GFLRB".
!         Then, we compute the first NCOFRB derivatives of GFLRB with
!         respect to a parameter (for us, REYNLD).  The first derivative
!         is stored in column 1 of RB, and so on.  Then we orthogonalize
!         the columns of RB.
!
!         We intend that NEQNFL >> NCOFRB, and RB is a matrix with orthogonal
!         columns, so that:
!
!           Transpose(RB) * RB = Identity(NCOFRB)
!
!
!         If GFL is any set of finite element coefficients, the corresponding
!         set of reduced basis coefficients can be computed as:
!
!           GRB = Transpose(RB) * (GFL-GFLRB)
!
!         If GRB is a set of reduced basis coefficients, a corresponding
!         set of finite element coefficients can be computed as:
!
!           GFL = GFLRB + RB * GRB.
!
!         While it is the case that you can expand and then reduce,
!         and always get the same result, it is not the case that
!         when you reduce and then expand you get the same result!
!
!         It is true, for ANY GRB, that
!
!           GRB = Transpose(RB) * RB * GRB
!
!         which follows from Transpose(RB) * RB = Identity(NCOFRB).
!
!         However, for a general GFL, it is the case that
!
!           GFL  = /= GFLRB + RB * Transpose(RB) * (GFL-GFLRB).
!
!         Only if GFL was generated from a reduced basis coefficient
!         vector will equality apply.  In other words, if GFL was generated
!         from a reduced basis coefficient:
!
!           GFL = GFLRB + RB * GRB
!
!         then
!
!           GFLRB + RB * Transpose(RB) * (GFL-GFLRB)
!           = GFLRB + RB * Transpose(RB) * (RB * GRB)
!           = GFLRB + RB *                       GRB
!           = GFL
!
!         so in this strictly limited case,
!
!           RB * Transpose(RB) = Identity(NEQNFL).
!
!    real ( kind = 8 ) RESRB(NCOFRB).
!    RESRB contains the residual in the reduced basis equations,
!    for the parameter values PAR and reduced basis coefficients GRB.
!
!    real ( kind = 8 ) SPLBMP(NPARB+2).
!    SPLBMP contains the spline coefficients for the bump.
!
!    real ( kind = 8 ) TAUBMP(NPARB+2).
!    TAUBMP contains the location of the spline abscissas for
!    the bump.  There are NPARB+2 of them, because the end values
!    of the spline are constrained to have particular values.
!
!    real ( kind = 8 ) TOLNEW.
!    TOLNEW is the convergence tolerance for the Newton
!    iteration.
!
!    real ( kind = 8 ) TOLOPT.
!    TOLOPT is the convergence tolerance for the optimization.
!    If TOLOPT is zero, then default values are used.
!
!    real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the Picard iteration.
!
!    real ( kind = 8 ) VOPT(LV).
!    VOPT provides real workspace for the optimization routines.
!
!    real ( kind = 8 ) WATEB.
!    WATEB is the multiplier of the bump control cost used
!    when computing the total cost.
!
!    real ( kind = 8 ) WATEP, WATEU, WATEV.
!    WATEP, WATEU and WATEV are weights used in computing the
!    cost function based on the costs of the flow discrepancy.
!
!    real ( kind = 8 ) XBL.
!    XBL is the X coordinate of the left corner of the bump.
!
!    real ( kind = 8 ) XBR.
!    XBR is the X coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) XOPT(MAXPAR).
!    XOPT is used by the optimization routines to hold only
!    the values of parameters which are allowed to vary.
!
!    real ( kind = 8 ) YBL.
!    YBL is the Y coordinate of the left corner of the bump.
!
!    real ( kind = 8 ) YBR.
!    YBR is the Y coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
  implicit none
!
!  Set parameters that are independent.
!
  integer ( kind = 4 ) liv
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm   
  integer ( kind = 4 ) maxnfl  
  integer ( kind = 4 ) maxnp              
  integer ( kind = 4 ) maxny
  integer ( kind = 4 ) maxpar     
  integer ( kind = 4 ) maxparb
  integer ( kind = 4 ) ncofrb
!
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,maxelm)
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) dopt(maxpar)
  real ( kind = 8 ) gflrb(maxnfl)
  real ( kind = 8 ) gfltar(maxnfl)
  real ( kind = 8 ) gfltmp(maxnfl)
  real ( kind = 8 ) grb(ncofrb)
  real ( kind = 8 ) grbarb(ncofrb)
  real ( kind = 8 ) grbopt(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(3,maxnp)
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) ipivrb(maxcofrb)
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivopt(liv)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) maxsim
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nopt
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nprof(2*maxny-1)
  integer ( kind = 4 ) numdif
  integer ( kind = 4 ) numopt
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(maxpar)
  real ( kind = 8 ) pararb(maxpar)
  real ( kind = 8 ) paropt(maxpar)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) resrb(maxcofrb)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) splbmp(maxparb+2)
  real ( kind = 8 ) taubmp(maxparb+2)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) vopt(lv)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xopt(maxpar)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(maxnp)
  real ( kind = 8 ) yrange
!
  ierror = 0
!
!  Copy the initial solution estimate.
!
  paropt(1:npar) = par(1:npar)
  grbopt(1:ncofrb) = grb(1:ncofrb)
!
!  Initialize the local optimization data.
!
  cost = 0.0D+00
  dopt(1:npar) = 1.0D+00
  ivopt(1:liv) = 0
  nopt = 0
  vopt(1:lv) = 0.0D+00
  xopt(1:maxpar) = 0.0D+00
!
!  Set the 611 data to default values,
!  and then modify some values.
!
  ival = 2
  call deflt(ival,ivopt,liv,lv,vopt)

  vopt(31) = tolopt
  vopt(32) = tolopt
  vopt(33) = tolopt
  vopt(34) = tolopt
  vopt(37) = tolopt

  ivopt(1) = 12
  ivopt(19) = 0
!
!  Set the step counters.
!
  numdif = 0
  numopt = 0
!
!  Take the next optimization step.
!
10    continue

  if ( maxopt < numopt ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'OPTDIFRB - Warning!'
    write ( *, '(a)' ) '  The number of optimization steps was exceeded.'
    return
  end if
!
!  Make a "copy" of PAR that only contains the free variables.
!
  nopt = 0
  do i = 1,npar
    if ( iopt(i) == 1) then
      nopt = nopt+1
      xopt(nopt) = paropt(i)
    end if
  end do
!
!  Call the optimizer to get a new parameter estimate.
!
  call snoit(dopt,cost,ivopt,liv,lv,nopt,vopt,xopt)
!
!  Copy the new free variable values back into PAR.
!
  nopt = 0
  do i = 1,npar
    if ( iopt(i) == 1) then
      nopt = nopt+1
      paropt(i) = xopt(nopt)
    end if
  end do
!
!  Apply Picard's method to the approximate solution GRBOPT.
!
    call picrb(arb,area,grbopt,ierror,indx,ipivrb,iwrite, &
      maxcofrb,maxelm,maxnfl,maxsim,nbcrb,ncofrb,nelem, &
      nferb,node,np,npar,nparf,nx,ny,paropt,phirb,rb,resrb,rmax, &
      tauflo,tolsim,xc,xrange,yc,yrange)
!
!  Apply Newton's method to the approximate solution GRBOPT.
!
    if ( tolnew < rmax ) then
      call newtrb(arb,area,grbopt,grbarb,ierror,indx,ipivrb, &
        iwrite,maxcofrb,maxelm,maxnew,maxnfl,nbcrb,ncofrb,nelem, &
        nferb,node,np,npar,nparf,nx,ny,paropt,pararb,phirb, &
        rb,resrb,rmax,tauflo,tolnew,xc,xrange,yc,yrange)

      if ( ierror /= 0) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'OPTDIFRB - Fatal error!'
        write ( *, '(a)' ) '  NEWTRB failed!'
        write ( *, '(a)' ) '  The parameters at which failure occurred:'
        write ( *, '(a)' ) ' '
        call prpar(iopt,npar,nparb,nparf,paropt)
        ierror = 1
        return
      end if
    end if
!
!  Compute the equivalent full basis solution GFLTMP = RB*GRB.
!
  call grb2fl ( gfltmp, gflrb, grbopt, maxnfl, ncofrb, neqnfl, rb )
!
!  Compute the cost function COST.
!
  call getcst(cost,costb,costp,costu,costv,gfltmp,gfltar, &
    indx,neqnfl,np,nparb,nprof,ny,splbmp, &
    taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)

  if ( ivopt(1) == 1 .and. 0 <= iwrite ) then
    call prpar(iopt,npar,nparb,nparf,paropt)
    write ( *, '(a,g14.6)' ) '  Cost = ',cost
  end if
!
!  If IVOPT(1) is 1, then this was a call for a legitimate
!  solution candidate.
!
!  If IVOPT(1) is 2, then this was a call for a temporary
!  solution used only for estimating the gradient.
!
!  Other values of IVOPT indicate acceptance or rejection
!  of the iteration.
!
  if ( ivopt(1) == 1) then
    numopt = numopt+1
  else if ( ivopt(1) == 2) then
    numdif = numdif+1
  else if ( 3 <= ivopt(1) .and. ivopt(1) <= 8 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Convergence to a minimizer was achieved!'
    return
  else if ( 8 < ivopt(1) ) then
    write ( *, '(a)' ) ' '
    write ( *, * ) '  Bad value of IVOPT(1) = ',ivopt(1)
    return
  end if

  go to 10

end
subroutine picfl(afl,area,eqn,gfl,ierror,ifs,indx,ipivfl,iwrite,ldafl, &
  maxsim,nelem,neqnfl,nlband,node,np,npar,nparf,numsim,par,phifl,region, &
  resfl,rmax,splflo,tauflo,tolsim,xc,xrange,yc,yrange)
!
!*****************************************************************************80
!
!! PICFL carries out simple iteration on the full Navier Stokes equations.
!
!
!  Discussion:
!
!    The simple iteration equations have the form:
!
!      Integral
!
!        dU/dx * dW/dx
!      + dU/dy * dW/dy
!      + reynld * (UOLD*dU/dx + VOLD*dU/dy + dP/dx) * W dx dy = 0
!
!      Integral
!
!        dV/dx * dW/dx
!      + dV/dy * dW/dy
!      + reynld * (UOLD*dV/dx + VOLD*dV/dy + dP/dy) * W dx dy = 0
!
!      Integral
!
!        (dU/dx + dV/dy) * Q dx dy = 0
!
!    Here W is a basis function for U and V, and Q is a basis
!    function for P.  UOLD and VOLD are the values of U and V
!    on the previous step of the iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) AFL(LDAFL,MAXNFL).
!    If Newton iteration is being carried out, AFL contains the
!    Jacobian matrix for the full system.
!    If Picard iteration is being carried out, AFL contains the
!    Picard matrix for the full system.
!
!    AFL is stored in LINPACK general band storage mode, with
!    logical dimensions (3*NLBAND+1, NEQNFL).
!
!    Where is the (I,J) entry of AFL actually stored?
!    AFL has actual storage for such an entry only if
!      -NLBAND <= I-J <= NLBAND.
!    In such a case, the (I,J) entry is actually stored in
!      AFL(I-J+2*NLBAND+1,J)
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!
!    or, if the element is isoperimetric,
!
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    character ( len = 2 ) EQN(MAXNFL).
!    EQN records the "type" of each equation that will be generated, and
!    which is associated with an unknown. 
!
!    'U'  A horizontal momentum equation.
!    'UB' The condition U = 0 applied at a node on the bump.
!    'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U = 0 applied at a node on a fixed wall.
!    'U0' A dummy value of U = 0 should be set.
!
!    'V'  A vertical momentum equation.
!    'VB' The condition V = 0 applied at a node on the bump.
!    'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V = 0 applied at a node on a fixed wall.
!    'V0' A dummy value of V = 0 should be set.
!
!    'P'  A continuity equation.
!    'PB' The condition P = 0 applied at (XMAX,YMAX).
!    'P0' A dummy value of P = 0 should be set.
!
!    real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
!
!    integer ( kind = 4 ) IERROR.
!    IERROR is an error flag.
!    0, no error occurred in this routine.
!    nonzero, an error occurred.
!
!    integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    Workspace, integer IPIVFL(NEQNFL).
!    IPIVFL is a pivot vector for the solution of the full
!    linear system.
!
!    integer ( kind = 4 ) IWRITE.
!    IWRITE controls the amount of output printed.
!    0, print out the least amount.
!    1, print out some.
!    2, print out a lot.
!
!    integer ( kind = 4 ) LDAFL.
!    LDAFL is the first dimension of the matrix AFL as declared in
!    the main program.  LDAFL must be at least 3*NLBAND+1.
!
!    integer ( kind = 4 ) MAXSIM.
!    MAXSIM is the maximum number of steps to take in one Picard
!    iteration.  A typical value is 20.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NLBAND.
!    NLBAND is the lower bandwidth of the matrix AFL.
!    The zero structure of AFL is assumed to be symmetric, and so
!    NLBAND is also the upper bandwidth of AFL.
! 
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!  
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
! 
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    integer ( kind = 4 ) NUMSIM.
!    NUMSIM is the number of simple iterations taken on a particular
!    call to the simple iteration routine.
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    real ( kind = 8 ) PHIFL(3,6,10,NELEM).
!    PHIFL contains the value of a finite element basis function, its
!    derivative, or other information, evaluated at the quadrature
!    points (which are the element midside nodes).
!
!    The meaning of the entry PHIFL(I,J,K,L) is as follows.
!    For the quadrature point I, and basis function J, in element L,
!    PHIFL(I,J,K,L) represents the value of:
!
!      K =  1, W, the finite element basis function for velocities;
!      K =  2, dWdX, the X derivative of W;
!      K =  3, dWdY, the Y derivative of W;
!      K =  4, Q, the finite element basis function for pressures;
!      K =  5, dQdX, the X derivative of Q;
!      K =  6, dQdY, the Y derivative of Q;
!      K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!      K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!      K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!      K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!    In particular, PHIFL(I,J,K,L) is the value of the quadratic
!    basis function W associated with local node J in element L,
!    evaluated at quadrature point I.
!
!    Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!    since there are only three linear basis functions.
!
!    character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottome, with tangential velocity specifications
!    there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    real ( kind = 8 ) RESFL(NEQNFL).
!    RESFL contains the residual in the full basis equations.
!
!    Output, real ( kind = 8 ) RMAX.
!    RMAX is the maximum absolute value of the entries of the residual
!    vector evaluated at the returned solution estimate.
!
!    real ( kind = 8 ) SPLFLO(NPARF).
!    SPLFLO contains the spline coefficients for the inflow.
!
!    real ( kind = 8 ) TAUFLO(NPARF).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow.  
!
!    real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the Picard iteration.
!
!    real ( kind = 8 ) XC(NP).
!    XC contains the X coordinates of the nodes.
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) afl(ldafl,neqnfl)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dxmax
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivfl(neqnfl)
  integer ( kind = 4 ) irmax
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) maxsim
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) numsim
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmax0
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax0
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
  ierror = 0
!
!  Get XMAX0, the norm of the initial guess GFL.
!
  ixmax = idamax ( neqnfl, gfl, 1 )
  xmax = abs(gfl(ixmax))
  xmax0 = xmax
!
!  Get RMAX0, the norm of the error RESFL of the initial guess, GFL.
!
  call fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar, &
    nparf,par,phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)

  irmax = idamax(neqnfl,resfl,1)
  rmax = abs(resfl(irmax))
  rmax0 = rmax

  numsim = 0

  if ( 2 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' Step	MxNorm(X)   IXmax  MxNorm(FX)	IRmax'
    write(*,'(i6,g14.6,i6,g14.6,i6)')numsim,xmax,ixmax,rmax,irmax
  end if

  if ( rmax0 < tolsim) then
    return
  end if
!
!  Do up to MAXSIM steps of simple iteration.
!
  do numsim = 1,maxsim
!
!  Get the simple iteration system matrix AFL evaluated at GFL.
!
    call picmfl(afl,area,eqn,gfl,indx,ldafl,nelem,neqnfl,nlband, &
      node,np,npar,par,phifl)
!
!  Factor the matrix.
!
    call dfacfl(afl,ldafl,neqnfl,nlband,nlband,ipivfl,info)

    if ( info /= 0) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PicFL - Fatal error!'
      write ( *, '(a)' ) '  The Picard matrix AFL is singular!'
      write ( *, * ) '  DFACFL returns nonzero INFO = ',info
      ierror = 1
      return
    end if
!
!  Get the right hand side, RESFL.
!
    call picvfl(eqn,ifs,indx,neqnfl,np,npar,nparf,par, &
      region,resfl,splflo,tauflo,xc,xrange,yc,yrange)

    ixmax = idamax(neqnfl,resfl,1)
    dxmax = abs(resfl(ixmax))
!
!  Solve the linear system AFL*GFL = RESFL.
!
    call dsolfl(afl,ldafl,neqnfl,nlband,nlband,ipivfl,resfl)
!
!  Compare RESFL and the previous estimate GFL.
!
    dxmax = 0.0D+00
    ixmax = 0
    do i = 1,neqnfl
      if ( dxmax <= abs ( resfl(i) - gfl(i) ) ) then
        ixmax = i
        dxmax = abs ( resfl(i) - gfl(i) )
      end if
    end do
!
!  Update GFL with the new estimate, and save its norm.
!
    do i = 1,neqnfl
      gfl(i) = resfl(i)
    end do

    ixmax = idamax(neqnfl,gfl,1)
    xmax = abs(gfl(ixmax))
!
!  Compute FX(GFL).
!
    call fxfl(area,eqn,gfl,ifs,indx,nelem,neqnfl,node,np,npar, &
      nparf,par,phifl,region,resfl,splflo,tauflo,xrange,yc,yrange)

    irmax = idamax(neqnfl,resfl,1)
    rmax = abs(resfl(irmax))
!
!  Print out
!
    if ( 2 <= iwrite ) then
      write(*,'(i6,g14.6,i6,g14.6,i6)')numsim,xmax,ixmax,rmax,irmax
    end if
!
!  Converged, Failed, or Continue?
!
    if ( rmax < tolsim * (rmax0+1.0D+00 )) then
      if ( 2 <= iwrite ) then
        write ( *, '(a)' ) 'PicFL - Residual acceptance.'
      end if
      return
    end if

    if ( 1000.0D+00 *rmax0 < rmax ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PicFL - Fatal error!'
      write ( *, * ) '  Simple iteration diverging on step ',numsim
      write ( *, * ) '  MxNorm of first X = ',xmax0
      write ( *, * ) '  MxNorm of last X =  ',xmax
      write ( *, * ) '  MxNorm of first FX  = ',rmax0
      write ( *, * ) '  MxNorm of last FX = ',rmax
      return
    end if

  end do

  if ( 3 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PicFL - Warning:'
    write ( *, '(a)' ) '  Simple iteration did not converge.'
    write ( *, * ) '  MxNorm of first X = ',xmax0
    write ( *, * ) '  MxNorm of last X =  ',xmax
    write ( *, * ) '  MxNorm of first FX  = ',rmax0
    write ( *, * ) '  MxNorm of last FX = ',rmax
  end if

  return
end
subroutine picmferb(arb,area,grb,maxcofrb,maxelm,nbcrb,ncofrb, &
  nelem,nferb,phirb,reynld)
!
!*****************************************************************************80
!
!! PICMFERB evaluates the simple iteration matrix for a reduced problem.
!
!
!  Discussion:
!
!    PICMFERB is given
!
!      GRB, the reduced basis coefficients of an approximate solution,
!      PHIRB, the reduced basis functions, evaluated at the quadrature
!        points,
!      REYNLD, the current Reynolds number,
!
!    and computes
!
!      ARB, the simple iteration matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) ARB(MAXCOFRB,MAXCOFRB).
!    ARB contains the Jacobian or Picard matrix for the reduced
!    Navier Stokes system, stored as an NCOFRB by NCOFRB array.
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXELM.
!    MAXELM is the maximum number of elements.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NFERB.
!    NFERB is the number of reduced basis coefficients that will
!    be determined via the finite element method.
!
!    real ( kind = 8 ) PHIRB(3,MAXCOFRB,15,MAXELM).
!    PHIRB contains the values of a finite element basis function
!    or its X or Y derivative, in a given element, at a given
!    quadrature point, for a particular reduced basis function.
!
!    For PHIRB(I,J,K,L), index J refers to the reduced basis
!    basis functions, for J = 0 to NCOFRB.
!
!    The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!      For the quadrature point I, and reduced basis function J,
!      in element L, PHIRB(I,J,K,L) represents the value of:
!
!        K = 1, WUrb, the finite element U velocity basis function;
!        K = 2, dWUrbdX, the X derivative of WUrb;
!        K = 3, dWUrbdY, the Y derivative of WUrb;
!        K = 4, WVrb, the finite element V velocity basis function;
!        K = 5, dWVrbdX, the X derivative of WVrb;
!        K = 6, dWVrbdY, the Y derivative of WVrb;
!        K = 7, Q, the finite element pressure basis function.
!        K = 8, dQrbdX, the X derivative of Qrb;
!        K = 9, dQrbdY, the Y derivative of Qrb.
!        K = 10, WU0rb, same as WUrb, with zero BC.
!        K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!        K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!        K = 13, WV0rb, same as WVrb, with zero BC.
!        K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!        K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the current value of the Reynolds number.
!    Normally, REYNLD is stored as PARA(NPARF+NPARB+1).
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
!
  real ( kind = 8 ) ar
  real ( kind = 8 ) arb(maxcofrb,maxcofrb)
  real ( kind = 8 ) area(3,maxelm)
  real ( kind = 8 ) dqjdx
  real ( kind = 8 ) dqjdy
  real ( kind = 8 ) dprbdx
  real ( kind = 8 ) dprbdy
  real ( kind = 8 ) durbdx
  real ( kind = 8 ) durbdy
  real ( kind = 8 ) dvrbdx
  real ( kind = 8 ) dvrbdy
  real ( kind = 8 ) dwu0dx
  real ( kind = 8 ) dwujdx
  real ( kind = 8 ) dwu0dy
  real ( kind = 8 ) dwujdy
  real ( kind = 8 ) dwv0dx
  real ( kind = 8 ) dwvjdx
  real ( kind = 8 ) dwv0dy
  real ( kind = 8 ) dwvjdy
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) icofrb
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) jcofrb
  logical s_eqi
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  real ( kind = 8 ) prb
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) urb
  real ( kind = 8 ) vrb
  real ( kind = 8 ) wu0
  real ( kind = 8 ) wv0
!
!  Zero out the FE rows of the matrix.
!
  do icofrb = nbcrb+1,nbcrb+nferb
    arb(icofrb,1:ncofrb) = 0.0D+00
  end do
!
!  Consider an element IELEM...
!
  do ielem = 1,nelem
!
!  ...and a quadrature point IQUAD...
!
    do iquad = 1,3

      ar = area(iquad,ielem)
!
!  For the given reduced coefficients GRB, and basis functions
!  PHIRB, evaluate U, V, and P, and their spatial derivatives.
!
      call uvpqrb(dprbdx,dprbdy,durbdx,durbdy,dvrbdx,dvrbdy,grb, &
        ielem,iquad,maxcofrb,maxelm,ncofrb,phirb,prb,urb,vrb)
!
!  Consider FE reduced basis function ICOFRB.
!
      do icofrb = nbcrb+1,nbcrb+nferb

        wu0    = phirb(iquad,icofrb,10,ielem)
        dwu0dx = phirb(iquad,icofrb,11,ielem)
        dwu0dy = phirb(iquad,icofrb,12,ielem)
        wv0    = phirb(iquad,icofrb,13,ielem)
        dwv0dx = phirb(iquad,icofrb,14,ielem)
        dwv0dy = phirb(iquad,icofrb,15,ielem)
!
!  Take the derivative with respect to basis function JCOFRB.
!
        do jcofrb = 1,ncofrb

          dwujdx = phirb(iquad,jcofrb,2,ielem)
          dwujdy = phirb(iquad,jcofrb,3,ielem)

          dwvjdx = phirb(iquad,jcofrb,5,ielem)
          dwvjdy = phirb(iquad,jcofrb,6,ielem)

          dqjdx  = phirb(iquad,jcofrb,8,ielem)
          dqjdy  = phirb(iquad,jcofrb,9,ielem)
!
!  The horizontal momentum equations.
!
          arb(icofrb,jcofrb) = arb(icofrb,jcofrb) &
               +ar*(dwujdx*dwu0dx + dwujdy*dwu0dy &
               +reynld*(urb*dwujdx+vrb*dwujdy+dqjdx)*wu0)
!
!  The vertical momentum equations.
!
          arb(icofrb,jcofrb) = arb(icofrb,jcofrb) &
               +ar*(dwvjdx*dwv0dx + dwvjdy*dwv0dy &
               +reynld*(urb*dwvjdx+vrb*dwvjdy+dqjdy)*wv0)

        end do
      end do
    end do
  end do

  return
end
subroutine picmfl ( afl, area, eqn, gfl, indx, ldafl, nelem, neqnfl, &
  nlband, node, np, npar, par, phifl )
!
!*****************************************************************************80
!
!! PICMFL computes the Picard iteration matrix for the full Navier Stokes equations.
!
!
!  The coefficients are:
!
!
!  d U-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + reynld * (Uold*dWj/dx+ Vold*dWj/dy) * Wi dx dy
!
!  d U-Eqn/d P-Coef:
!
!    Integral
!
!    reynld * dQj/dx * Wi dx dy
!
!  d V-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dx * dWi/dx + dWj/dy * dWi/dy
!    + reynld * (Uold*dWj/dx + Vold*dWj/dy) * Wi dx dy
!
!  d V-Eqn/d P-Coef:
!
!    Integral
!
!    reynld * dQj/dy * Wi dx dy
!
!  d P-Eqn/d U-Coef:
!
!    Integral
!
!      dWj/dx * Qi dx dy
!
!  d P-Eqn/d V-Coef:
!
!    Integral
!
!      dWj/dy * Qi dx dy
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) A(LDAFL,NEQNFL), contains the
!    coefficients of the Picard iteration matrix.
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!  EQN    Input, character ( len = 2 ) EQN(NEQNFL).
!         EQN records the "type" of each equation that will be generated, and
!         which is associated with an unknown.  Note that most boundary
!         conditions do not result in an equation.  The current values are:
!
!         'U'  The horizontal momentum equation.
!         'UB' The condition U = 0 applied at a node on the bump.
!         'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!         'UW' The condition U = 0 applied at a node on a fixed wall.
!
!         'V'  The vertical momentum equation.
!         'VB' The condition V = 0 applied at a node on the bump.
!         'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!         'VW' The condition V = 0 applied at a node on a fixed wall.
!
!         'P'  The continuity equation.
!         'PB' The condition P = 0 applied at (XMAX,YMAX).
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!         G is the current solution vector, in which are stored
!         the finite element coefficients that define the velocity
!         and pressure functions, U, V and P.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!    Input, integer ( kind = 4 ) LDAFL, the first dimension of the matrix AFL.
!
!    Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!    Input, integer ( kind = 4 ) NEQNFL, the number of finite element equations used
!    to define the horizontal and vertical velocities and the
!    pressure.
!
!    Input, integer ( kind = 4 ) NLBAND.
!    The lower bandwidth of the matrix A.  The zero structure of A
!    is assumed to be symmetric, and so NLBAND is also the upper
!    bandwidth of A. 
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM).
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  NPAR   Input, integer ( kind = 4 ) NPAR.
!         The number of parameters.  NPAR = NPARF + NPARB + 1.
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!         The number of parameters associated with the position and
!         shape of the bump.
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1. 
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!         PAR is the current set of parameter values, including the
!         Reynolds parameter, the flow parameters, and the bump parameters.
!
!  PHIFL  Input, real ( kind = 8 ) PHIFL(3,6,10,NELEM). 
!         PHIFL contains the value of a finite element basis function, its
!         derivative, or other information, evaluated at the quadrature
!         points.
!         The meaning of the entry PHIFL(I,J,K,L) is as follows. 
!         For the quadrature point I, and basis function J, in element L,
!         PHIFL(I,J,K,L) represents the value of:
!
!           K =  1, W, the finite element basis function for velocities;
!           K =  2, dWdX, the X derivative of W;
!           K =  3, dWdY, the Y derivative of W;
!           K =  4, Q, the finite element basis function for pressures;
!           K =  5, dQdX, the X derivative of Q;
!           K =  6, dQdY, the Y derivative of Q;
!           K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!           K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!           K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!           K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!         In particular, PHIFL(I,J,K,L) is the value of the quadratic
!         basis function W associated with local node J in element L,
!         evaluated at quadrature point I.
!
!         Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!         since there are only three linear basis functions.
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
!
  real ( kind = 8 ) afl(ldafl,neqnfl)
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dqjdx
  real ( kind = 8 ) dqjdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwidx
  real ( kind = 8 ) dwidy
  real ( kind = 8 ) dwjdx
  real ( kind = 8 ) dwjdy
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iuse
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhor
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) jprs
  integer ( kind = 4 ) jq
  integer ( kind = 4 ) jver
  logical s_eqi
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) qi
  real ( kind = 8 ) reynld
  real ( kind = 8 ) term
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) wi
!
  reynld = par(npar)

  do i = 1,3*nlband+1
    afl(i,1:neqnfl) = 0.0D+00
  end do
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1,nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1,3

      ar = area(iquad,ielem)
!
!  Evaluate U, V and P at the IQUAD-th quadrature point.
!
      call uvpqfl(dpdx,dpdy,dudx,dudy,dvdx,dvdy,gfl,ielem,indx, &
        iquad,nelem,neqnfl,node,np,p,phifl,u,v)
!
!  Consider each node in the element.
!
      do iq = 1,6

        ip = node(iq,ielem)

        wi = phifl(iquad,iq,1,ielem)
        dwidx = phifl(iquad,iq,2,ielem)
        dwidy = phifl(iquad,iq,3,ielem)
        qi = phifl(iquad,iq,4,ielem)

        ihor = indx(1,ip)
        iver = indx(2,ip)
        iprs = indx(3,ip)
!
!  Now compute the derivatives of the functions associated
!  with U, V and P, with respect to the coefficients associated
!  with basis vectors at each node of the element.
!
        do jq = 1,6

          jp = node(jq,ielem)

          dwjdx = phifl(iquad,jq,2,ielem)
          dwjdy = phifl(iquad,jq,3,ielem)

          dqjdx = phifl(iquad,jq,5,ielem)
          dqjdy = phifl(iquad,jq,6,ielem)

          jhor = indx(1,jp)
          jver = indx(2,jp)
          jprs = indx(3,jp)
!
!  Contributions of the JHOR horizontal velocity to the U, V, and
!  P equations.
!
          if ( eqn(ihor) == 'U') then

            term = ar*(dwjdx*dwidx+dwjdy*dwidy+ &
                 reynld*(u*dwjdx+v*dwjdy)*wi)

            iuse = ihor-jhor+2*nlband+1
            afl(iuse,jhor) = afl(iuse,jhor)+term

          end if

          if ( 0 < iprs ) then
            if ( eqn(iprs) == 'P') then
              term = ar*dwjdx*qi
              iuse = iprs-jhor+2*nlband+1
              afl(iuse,jhor) = afl(iuse,jhor)+term
            end if
          end if
!
!  Contributions of the JVER vertical velocity variable to the
!  U, V and P equations.
!
          if ( eqn(iver) == 'V') then

            term = ar*(dwjdx*dwidx+dwjdy*dwidy &
                 +reynld*(u*dwjdx+v*dwjdy)*wi)

            iuse = iver-jver+2*nlband+1
            afl(iuse,jver) = afl(iuse,jver)+term
          end if

          if ( 0 < iprs ) then
            if ( eqn(iprs) == 'P') then
              term = ar*dwjdy*qi
              iuse = iprs-jver+2*nlband+1
              afl(iuse,jver) = afl(iuse,jver)+term
            end if
          end if
!
!  Contributions of the JPRS pressure to the U and V equations.
!
          if ( 0 < jprs ) then

            if ( eqn(ihor) == 'U') then
              term = ar*reynld*dqjdx*wi
              iuse = ihor-jprs+2*nlband+1
              afl(iuse,jprs) = afl(iuse,jprs)+term
            end if

            if ( eqn(iver) == 'V') then
              term = ar*reynld*dqjdy*wi
              iuse = iver-jprs+2*nlband+1
              afl(iuse,jprs) = afl(iuse,jprs)+term
            end if

          end if

        end do
      end do
    end do
  end do
!
!  Set up the equations that enforce boundary conditions.
!
  do ip = 1,np

    ihor = indx(1,ip)
    iver = indx(2,ip)
    iprs = indx(3,ip)

    if ( eqn(ihor) == 'UB'.or.eqn(ihor) == 'UI'.or. &
         eqn(ihor) == 'UW'.or. eqn(ihor) == 'U0') then
      afl(2*nlband+1,ihor) = 1.0D+00
    end if

    if ( eqn(iver) == 'VB'.or. eqn(iver) == 'VI'.or. &
         eqn(iver) == 'VW'.or. eqn(iver) == 'V0') then
      afl(2*nlband+1,iver) = 1.0D+00
    end if

    if ( 0 < iprs ) then
      if ( eqn(iprs) == 'PB') then
        afl(2*nlband+1,iprs) = 1.0D+00
      else if ( eqn(iprs) == 'P0') then
        afl(2*nlband+1,iprs) = 1.0D+00
      end if
    end if

  end do

  return
end
subroutine picrb ( arb, area, grb, ierror, indx, ipivrb, iwrite, &
  maxcofrb, maxelm, maxnfl, maxsim, nbcrb, ncofrb, nelem, nferb, node, &
  np, npar, nparf, nx, ny, par, phirb, rb, resrb, rmax, tauflo, &
  tolsim, xc, xrange, yc, yrange )
!
!*****************************************************************************80
!
!! PICRB carries out simple iteration on the reduced Navier Stokes equations.
!
!
!  Discussion:
!
!    The simple iteration equations have the form:
!
!    Integral
!
!      dU/dx * dW/dx + dU/dy * dW/dy
!    + reynld * (URB*dU/dx + VRB*dU/dy + dP/dx) * W dx dy = 0
!
!    Integral
!
!      dV/dx * dW/dx + dV/dy * dW/dy
!    + reynld * (URB*dV/dx + VRB*dV/dy + dP/dy) * W dx dy = 0
!
!    Here W is a basis function for U and V.
!    UOLD and VOLD are the values of U and V
!    on the previous step of the iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Workspace, real ( kind = 8 ) ARB(MAXNRB,MAXNRB).
!    ARB contains the Jacobian or Picard matrix for the reduced
!    Navier Stokes system, stored as a dense NCOFRB by NCOFRB array.
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!  GRB    Input, real ( kind = 8 ) GRB(NCOFRB).
!         GRB contains the reduced basis coefficients of the current
!         estimate of the state solution.
!
!  GRBTMP Workspace, real ( kind = 8 ) GRBTMP(NCOFRB).
!
!  IERROR Output, integer ( kind = 4 ) IERROR.
!         0, no error occurred.
!         1, an error occurred.  The matrix was singular.
!
!  IPIVRB Workspace, integer IPIVRB(NCOFRB).
!
!  IWRITE Input, integer ( kind = 4 ) IWRITE.
!         IWRITE controls the amount of output printed.
!
!  MAXSIM Input, integer ( kind = 4 ) MAXSIM.
!         MAXSIM is the maximum number of simple iteration steps
!         that may be taken.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB.
!         NCOFRB is the number of basis functions, reduced state equations and
!         coefficients in the reduced basis system.
!
!  NPAR   Input, integer ( kind = 4 ) NPAR.
!         The number of parameters.  NPAR = NPARF + NPARB + 1.
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!         The number of parameters associated with the position and
!         shape of the bump.
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1. 
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!         PAR is the current set of parameter values, including the
!         Reynolds parameter, the flow parameters, and the bump parameters.
!
!  PHIRB  Input, real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    Workspace, real ( kind = 8 ) RESRB(NCOFRB).
!
!    Input, real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the iteration.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) arb(maxcofrb,ncofrb)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dxmax
  real ( kind = 8 ) grb(ncofrb)
  real ( kind = 8 ) grbtmp(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivrb(ncofrb)
  integer ( kind = 4 ) irmax
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) ixmax
  integer ( kind = 4 ) maxsim
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) numsim
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) rb(maxnfl,maxcofrb)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) rmax
  real ( kind = 8 ) rmax0
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmax0
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
  if ( ncofrb <= 0) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PICRB - Fatal error!'
    write ( *, '(a,i6)' ) '  NCOFRB <= 0, NCOFRB=',ncofrb
    stop
  end if

  reynld = par(npar)

  ierror = 0
!
!  Get XMAX0, the norm of the initial guess GRB.
!
  ixmax = idamax(ncofrb,grb,1)
  xmax = abs(grb(ixmax))
  xmax0 = xmax
!
!  Get RMAX0, the norm of the error RESRB of the initial guess, GRB.
!
  call fxrb(area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
    nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb, &
    resrb,reynld,tauflo,xc,xrange,yc,yrange)

  irmax = idamax(ncofrb,resrb,1)
  rmax = abs(resrb(irmax))
  rmax0 = rmax

  numsim = 0

  if ( 2 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' Step   MxNorm(X)   IXmax  MxNorm(FX)   IRmax'
    write(*,'(i6,g14.6,i6,g14.6,i6)')numsim,xmax,ixmax,rmax,irmax
  end if

  if ( rmax0 < tolsim) then
    return
  end if
!
!  Do up to MAXSIM steps of simple iteration.
!
  do numsim = 1,maxsim
!
!  Get the simple iteration system matrix ARB evaluated at GRB.
!
    call fpbcrb(arb,indx,maxcofrb,maxnfl,nbcrb,ncofrb, &
      nelem,node,np,nx,ny,rb,xc,xrange,yc,yrange)

    call picmferb(arb,area,grb,maxcofrb,maxelm,nbcrb,ncofrb, &
      nelem,nferb,phirb,reynld)
!
!  Factor the matrix.
!
    call dfacrb(arb,maxcofrb,ncofrb,ipivrb,info)

    if ( info /= 0) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PicRB - Fatal error!'
      write ( *, '(a)' ) '  The Picard matrix ARB is singular!'
      write ( *, * ) '  DFACRB returns nonzero INFO = ',info
      ierror = 1
      return
    end if
!
!  Get the right hand side, RESRB.
!
    grbtmp(1:ncofrb) = 0.0D+00

    call fxbcrb(grbtmp,indx,maxcofrb,maxnfl,nbcrb,ncofrb,nelem, &
      node,np,npar,nparf,nx,ny,par,rb,resrb,tauflo,xc,xrange,yc,yrange)

    do i = 1,nbcrb
      resrb(i) = -resrb(i)
    end do

    call picvferb(area,grb,maxcofrb,maxelm,nbcrb, &
      ncofrb,nelem,nferb,npar,par,phirb,resrb)
!
!  Solve the linear system ARB*GRB = RESRB.
!
    call dsolrb(arb,maxcofrb,ncofrb,ipivrb,resrb)
!
!  Compare RESRB and the previous estimate GRB.
!
    dxmax = 0.0D+00
    ixmax = 0
    do i = 1,ncofrb
      if ( dxmax <= abs ( resrb(i) - grb(i) ) ) then
        ixmax = i
        dxmax = abs ( resrb(i) - grb(i) )
      end if
    end do
!
!  Update GRB with the new estimate, and save its norm.
!
    grb(1:ncofrb) = resrb(1:ncofrb)

    ixmax = idamax(ncofrb,grb,1)
    xmax = abs(grb(ixmax))
!
!  Compute FX(GRB).
!
    call fxrb(area,grb,indx,maxcofrb,maxelm,maxnfl,nbcrb,ncofrb, &
      nelem,nferb,node,np,npar,nparf,nx,ny,par,phirb,rb, &
      resrb,reynld,tauflo,xc,xrange,yc,yrange)

    irmax = idamax(ncofrb,resrb,1)
    rmax = abs(resrb(irmax))
!
!  Print out
!
    if ( 2 <= iwrite ) then
      write(*,'(i6,g14.6,i6,g14.6,i6)')numsim,xmax,ixmax,rmax,irmax
    end if
!
!  Converged, Failed, or Continue?
!
    if ( rmax < tolsim * (rmax0+1.0D+00 )) then
      if ( 2 <= iwrite ) then
        write ( *, '(a)' ) 'PicRB - Residual acceptance.'
      end if
      return
    end if

    if ( 1000.0D+00 *rmax0 < rmax ) then
      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PICRB - Fatal error!'
      write ( *, * ) '  Simple iteration diverging on step ',numsim
      write ( *, * ) '  MxNorm of first X = ',xmax0
      write ( *, * ) '  MxNorm of last X =  ',xmax
      write ( *, * ) '  MxNorm of first FX  = ',rmax0
      write ( *, * ) '  MxNorm of last FX = ',rmax
      return
    end if

  end do

  if ( 3 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PICRB - Warning:'
    write ( *, '(a)' ) '  Simple iteration did not converge.'
    write ( *, * ) '  MxNorm of first X = ',xmax0
    write ( *, * ) '  MxNorm of last X =  ',xmax
    write ( *, * ) '  MxNorm of first FX  = ',rmax0
    write ( *, * ) '  MxNorm of last FX = ',rmax
  end if

  return
end
subroutine picvferb(area,grb,maxcofrb,maxelm,nbcrb,ncofrb,nelem,nferb, &
  npar,par,phirb,resrb)
!
!*****************************************************************************80
!
!! PICVFERB computes the finite element portion of the right hand
!  side of the Picard iteration for the reduced Navier Stokes
!  equations.
!
!  Discussion:
!
!    The right hand side is simply the basis solution GFLRB
!    multiplied by the iteration matrix and negated. 
!    The easiest way to access the solution in GFLRB is to set
!    a temporary copy of GRB to zero, and call UVPQRB.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Input, real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!    Input, integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of basis functions, reduced state equations and
!    coefficients in the reduced basis system.
!
!    Input, integer ( kind = 4 ) NPAR.
!    The number of parameters.  NPAR = NPARF + NPARB + 1.
!    The parameters control the shape of the inflow,
!    the shape of the bump obstacle, and the strength of the
!    flow.
!
!    Input, integer ( kind = 4 ) NPARB.
!    The number of parameters associated with the position and
!    shape of the bump.
!    Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!    Input, integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1. 
!
!    Input, real ( kind = 8 ) PAR(NPAR).
!    PAR is the current set of parameter values, including the
!    Reynolds parameter, the flow parameters, and the bump parameters.
!
!  PHIRB  Input, real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!    Output, real ( kind = 8 ) RESRB(NCOFRB).
!    For this routine, RESRB returns the right hand side of
!    the Picard iteration system.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) npar
!
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dprbdx
  real ( kind = 8 ) dprbdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) durbdx
  real ( kind = 8 ) durbdy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dvrbdx
  real ( kind = 8 ) dvrbdy
  real ( kind = 8 ) dwuidx
  real ( kind = 8 ) dwuidy
  real ( kind = 8 ) dwvidx
  real ( kind = 8 ) dwvidy
  real ( kind = 8 ) grb(ncofrb)
  real ( kind = 8 ) grbtmp(ncofrb)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) icofrb
  integer ( kind = 4 ) iquad
  logical s_eqi
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) nferb
  real ( kind = 8 ) p
  real ( kind = 8 ) prb
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) resrb(ncofrb)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) u
  real ( kind = 8 ) urb
  real ( kind = 8 ) v
  real ( kind = 8 ) vrb
  real ( kind = 8 ) wui
  real ( kind = 8 ) wvi
!
  reynld = par(npar)
  grbtmp(1:ncofrb) = 0.0D+00

  do icofrb = nbcrb+1,nbcrb+nferb
    resrb(icofrb) = 0.0D+00
  end do
!
!  Approximate the integral by summing over all elements.
!
  do ielem = 1,nelem
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1,3

      ar = area(iquad,ielem)
!
!  Evaluate the full solution GFLRB at which the reduced basis
!  was generated.  This is the implicit "1" coefficient in the
!  set of reduced basis coefficients, which must be multiplied
!  by the Picard coefficients and carried to the right hand side.
!
!  We do this by using GRBTMP, set to 0.
!
      call uvpqrb(dpdx,dpdy,dudx,dudy,dvdx,dvdy,grbtmp, &
        ielem,iquad,maxcofrb,maxelm,ncofrb,phirb,p,u,v)
!
!  Evaluate the reduced basis solution GRB from the previous iterate.
!
      call uvpqrb(dprbdx,dprbdy,durbdx,durbdy,dvrbdx,dvrbdy,grb, &
        ielem,iquad,maxcofrb,maxelm,ncofrb,phirb,prb,urb,vrb)
!
!  Now consider each reduced basis function, and retrieve the
!  corresponding values of the U and V basis functions.
!
      do icofrb = nbcrb+1,nbcrb+nferb

        wui    = phirb(iquad,icofrb,1,ielem)
        dwuidx = phirb(iquad,icofrb,2,ielem)
        dwuidy = phirb(iquad,icofrb,3,ielem)

        wvi    = phirb(iquad,icofrb,4,ielem)
        dwvidx = phirb(iquad,icofrb,5,ielem)
        dwvidy = phirb(iquad,icofrb,6,ielem)
!
!  The horizontal velocity equations.
!
        resrb(icofrb) = resrb(icofrb) &
             -ar*(dudx*dwuidx + dudy*dwuidy &
             +reynld*(urb*dudx+vrb*dudy+dpdx)*wui)
!
!  The vertical velocity equations.
!
        resrb(icofrb) = resrb(icofrb) &
             -ar*(dvdx*dwvidx + dvdy*dwvidy &
             +reynld*(urb*dvdx+vrb*dvdy+dpdy)*wvi )

      end do
    end do
  end do

  return
end
subroutine picvfl(eqn,ifs,indx,neqnfl,np,npar,nparf,par, &
  region,resfl,splflo,tauflo,xc,xrange,yc,yrange)
!
!*****************************************************************************80
!
!! PICVFL computes the Picard right hand side for the full equations.
!
!
!  Discussion:
!
!    The Picard iteration equations have the form:
!
!    Integral
!
!      dU/dx * dW/dx
!    + dU/dy * dW/dy
!    + reynld * (UOLD*dU/dx + VOLD*dU/dy + dP/dx) * W dx dy = 0
!
!    Integral
!
!      dV/dx * dW/dx
!    + dV/dy * dW/dy
!    + reynld * (UOLD*dV/dx + VOLD*dV/dy + dP/dy) * W dx dy = 0
!
!    Integral
!
!      (dU/dx + dV/dy) * Q dx dy = 0
!
!    Here W is a basis function for U and V, and Q is a basis
!    function for P.  UOLD and VOLD are the values of U and V
!    on a previous step of the iteration.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  character ( len = 2 ) eqn(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iver
  real ( kind = 8 ) par(npar)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) ubc
  real ( kind = 8 ) vbc
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xval
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
  real ( kind = 8 ) yval
!
!  Initialize the right hand side to zero.
!
  do i = 1,neqnfl
    resfl(i) = 0.0D+00
  end do
!
  do ip = 1,np

    xval = xc(ip)
    yval = yc(ip)

    ihor = indx(1,ip)
    iver = indx(2,ip)

    if ( eqn(ihor) == 'UI'.or.eqn(iver).eq.'VI') then

      call flowbc(ifs,npar,nparf,par,region,splflo,tauflo, &
        ubc,vbc,xrange,xval,yrange,yval)

      if ( eqn(ihor) == 'UI') then
        resfl(ihor) = ubc
      end if

      if ( eqn(iver) == 'VI') then
        resfl(iver) = vbc
      end if

    end if

  end do

  return
end
subroutine reysen(area,eqn,indx,isen,maxcofrb,maxnfl,nelem, &
  neqnfl,node,np,phifl,resfl,reynld,senfl)
!
!*****************************************************************************80
!
!! REYSEN sets up a right hand side for a REYNLD sensitivity equation.
!
!
!  Discussion:
!
!    The routine sets up the right hand side vector associated with 
!    the ISEN-th order sensitivities with respect to the REYNLD parameter 
!    of a given state function (U,V,P).
!
!    In order to compute the right hand side for the ISEN-th order,
!    a state solution and the sensitivities for all orders less than
!    ISEN must already have been computed.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) AREA(3,MAXELM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    character ( len = 2 ) EQN(MAXNFL).
!    EQN records the "type" of each equation that will be generated, and
!    which is associated with an unknown. 
!
!    'U'  A horizontal momentum equation.
!    'UB' The condition U = 0 applied at a node on the bump.
!    'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U = 0 applied at a node on a fixed wall.
!    'U0' A dummy value of U = 0 should be set.
!    'V'  A vertical momentum equation.
!    'VB' The condition V = 0 applied at a node on the bump.
!    'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V = 0 applied at a node on a fixed wall.
!    'V0' A dummy value of V = 0 should be set.
!    'P'  A continuity equation.
!    'PB' The condition P = 0 applied at (XMAX,YMAX).
!    'P0' A dummy value of P = 0 should be set.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    Input, integer ( kind = 4 ) ISEN.
!    ISEN is the order of the sensitivity to be calculated.
!
!    integer ( kind = 4 ) MAXCOFRB.
!    MAXCOFRB is the maximum legal value for NCOFRB, the number
!    of coefficients used to specify a particular reduced basis
!    solution.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    real ( kind = 8 ) PHIFL(3,6,10,NELEM).
!    PHIFL contains the value of a finite element basis function, its
!    derivative, or other information, evaluated at the quadrature
!    points (which are the element midside nodes).
!    The meaning of the entry PHIFL(I,J,K,L) is as follows.
!    For the quadrature point I, and basis function J, in element L,
!    PHIFL(I,J,K,L) represents the value of:
!
!      K =  1, W, the finite element basis function for velocities;
!      K =  2, dWdX, the X derivative of W;
!      K =  3, dWdY, the Y derivative of W;
!      K =  4, Q, the finite element basis function for pressures;
!      K =  5, dQdX, the X derivative of Q;
!      K =  6, dQdY, the Y derivative of Q;
!      K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!      K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!      K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!      K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!    In particular, PHIFL(I,J,K,L) is the value of the quadratic
!    basis function W associated with local node J in element L,
!    evaluated at quadrature point I.
!    Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!    since there are only three linear basis functions.
!
!    real ( kind = 8 ) RESFL(NEQNFL).
!    RESFL contains the residual in the full basis equations.
!
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the value of the Reynolds number.
!
!    real ( kind = 8 ) SENFL(MAXNFL,MAXCOFRB).
!    Columns 1 through NSENFL of SENFL contain the sensitivities
!    of the full solution with respect to the REYNLD parameter, for
!    orders 0 through NSENFL-1.
!    SENFL(I,J) contains the (J-1)-th sensitivity of the I-th full unknown
!    with respect to REYNLD.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) ar
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) dpdx(maxcofrb)
  real ( kind = 8 ) dpdy(maxcofrb)
  real ( kind = 8 ) dudx(maxcofrb)
  real ( kind = 8 ) dudy(maxcofrb)
  real ( kind = 8 ) dvdx(maxcofrb)
  real ( kind = 8 ) dvdy(maxcofrb)
  character ( len = 2 ) eqn(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) isen
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) jsen
  integer ( kind = 4 ) jsendx
  integer ( kind = 4 ) nbinom
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p(maxcofrb)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) reynld
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) term
  real ( kind = 8 ) u(maxcofrb)
  real ( kind = 8 ) v(maxcofrb)
  real ( kind = 8 ) wi
!
!  Check the value of REYNLD.
!
  if ( reynld <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REYSEN - Fatal error!'
    write ( *, * ) '  Nonpositive value of REYNLD = ',reynld
    stop
  end if
!
!  Check the value of ISEN.
!
  if ( isen <= 0) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REYSEN - Fatal error!'
    write ( *, * ) '  The input value of ISEN is ',isen
    write ( *, '(a)' ) '  but ISEN must be strictly positive.'
    stop
  end if

  if ( maxcofrb-1 < isen ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'REYSEN - Fatal error!'
    write ( *, * ) '  The input value of ISEN is ',isen
    write ( *, * ) '  but the limit is MAXCOFRB-1 = ',maxcofrb-1
    stop
  end if
!
!  Zero out the right hand side vector RESFL.
!
  do i = 1,neqnfl
    resfl(i) = 0.0D+00
  end do
!
!  Approximate the integrated ISEN-th order sensitivity equations by
!  adding the contribution from element IELEM.
!
  do ielem = 1,nelem
!
!  In element IELEM, approximate the integrals by moving to
!  quadrature point IQUAD. 
!
    do iquad = 1,3

      ar = area(iquad,ielem)
!
!  We are about to compute the sensitivity of order ISEN, which
!  is stored in vector entry ISEN+1, or matrix column ISEN+1.
!  Evaluate the fundamental solution (U,V,P), and the first ISEN-1
!  sensitivities, (U',V',P'), (U'',V'',P''), and so on, stored in SENFL.
!
      do jsen = 0,isen-1
        jsendx = jsen+1
        call uvpqfl(dpdx(jsendx),dpdy(jsendx),dudx(jsendx), &
             dudy(jsendx),dvdx(jsendx),dvdy(jsendx), &
             senfl(1,jsendx),ielem,indx,iquad,nelem,neqnfl,node, &
             np,p(jsendx),phifl,u(jsendx),v(jsendx))
      end do
!
!  Now consider a node with local index IQ, and global index IP,
!  whose quadratic basis function evaluated at the quadrature point
!  IQUAD has value WI.  Evaluate the right hand sides of equations
!  IHOR and IVER and add the contributions to the total.
!
      do iq = 1,6

        ip = node(iq,ielem)
        wi = phifl(iquad,iq,1,ielem)
        ihor = indx(1,ip)
        iver = indx(2,ip)

        if ( eqn(ihor) == 'U') then

          term = 0.0D+00

          do jsen = 1,isen-1
            jsendx = jsen+1
            term = term+reynld*nbinom(isen,jsen)*(u(isen-jsen+1)*dudx(jsendx) &
              +v(isen-jsen+1)*dudy(jsendx))
          end do

          do jsen = 0,isen-1
            jsendx = jsen+1
            term = term+isen*nbinom(isen-1,jsen)*(u(isen-jsen)*dudx(jsendx) &
                 +v(isen-jsen)*dudy(jsendx))
          end do

          term = term+isen*dpdx(isen)

          resfl(ihor) = resfl(ihor)-ar*term*wi

        end if

!
!  Note that the vertical right hand side should be obtainable
!  from the horizontal right hand side by interchanging U and V,
!  and X and Y.
!
        if ( eqn(iver) == 'V') then

          term = 0.0D+00

          do jsen = 1,isen-1
            jsendx = jsen+1
            term = term+reynld*nbinom(isen,jsen)*(v(isen-jsen+1)*dvdy(jsendx) &
               +u(isen-jsen+1)*dvdx(jsendx))
          end do

          do jsen = 0,isen-1
            jsendx = jsen+1
            term = term+isen*nbinom(isen-1,jsen)*(v(isen-jsen)*dvdy(jsendx) &
                +u(isen-jsen)*dvdx(jsendx))
          end do

          term = term+isen*dpdy(isen)

          resfl(iver) = resfl(iver)-ar*term*wi

        end if

      end do
    end do
  end do

  return
end
subroutine test2(gfl,grb,ihi,ilo,indx,maxcofrb,maxelm,ncofrb, &
  nelem,neqnfl,node,np,phifl,phirb)
!
!*****************************************************************************80
!
!! TEST2 compares full and reduced U, V, and P at the quadrature points.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dprbdx
  real ( kind = 8 ) dprbdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) durbdx
  real ( kind = 8 ) durbdy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dvrbdx
  real ( kind = 8 ) dvrbdy
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi2
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo2
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p
  real ( kind = 8 ) prb
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) u
  real ( kind = 8 ) urb
  real ( kind = 8 ) v
  real ( kind = 8 ) vrb
!
  if ( ilo < 1) then
    ilo2 = 1
  else
    ilo2 = ilo
  end if

  if ( nelem < ihi ) then
    ihi2 = nelem
  else
    ihi2 = ihi
  end if
!
!  Consider an element IELEM...
!
  write ( *, '(a)' ) ' '
  write ( *, * ) 'Elements from ILO2 = ',ilo2,' to IHI2=',ihi2

  do ielem = ilo2,ihi2
!
!  ...and a quadrature point IQUAD...
!
    do iquad = 1,3
!
!  Evaluate U, V, and P for GFL and for GRB.
!
      call uvpqfl(dpdx,dpdy,dudx,dudy,dvdx,dvdy,gfl,ielem,indx, &
           iquad,nelem,neqnfl,node,np,p,phifl,u,v)

      call uvpqrb(dprbdx,dprbdy,durbdx,durbdy,dvrbdx,dvrbdy,grb, &
           ielem,iquad,maxcofrb,maxelm,ncofrb,phirb,prb,urb,vrb)

      write ( *, '(a)' ) ' '
      write ( *, * ) '  Element ',ielem
      write ( *, * ) '  Quad point ',iquad
      write(*,'(a,3g14.6)')'Full U,V,P:   ',u,v,p
      write(*,'(a,3g14.6)')'Reduced U,V,P:',urb,vrb,prb

    end do
  end do

  return
end
subroutine test3 ( maxcofrb, maxnfl, ncofrb, neqnfl, rb, senfl, senrb )
!
!*****************************************************************************80
!
!! TEST3 verifies that RB*RFACT = SENFL
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) neqnfl
!
  real ( kind = 8 ) dmax
  integer ( kind = 4 ) i
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) k
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) senrb(maxcofrb,maxcofrb)
  real ( kind = 8 ) temp
!
  dmax = 0.0D+00  
  imax = 0
  jmax = 0

  do i = 1,neqnfl
    do j = 1,ncofrb
      temp = 0.0D+00
      do k = 1,ncofrb
        temp = temp+rb(i,k)*senrb(k,j)
      end do
      temp = abs(temp-senfl(i,j))
      if ( dmax <= temp ) then
        dmax = temp
        imax = i
        jmax = j
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST3 - Results:'
  write ( *, '(a)' ) '  The maximum difference between RB*SENRB and SENFL'
  write ( *, * ) '  is ',dmax
  write ( *, * ) ' for (I,J) = ',imax,jmax

  return
end
subroutine test4(afl,area,difcof,dpar,drey,eqn,gfl,gflafl, &
  ifs,ijac,indx,ipar,ipivfl,iwrite,ldafl,maxcofrb,maxelm, &
  maxnew,maxnfl,ncofrb,nelem,neqnfl,nlband,node,np,npar, &
  nparf,nsenfl,par,parafl,phifl,region,resfl, &
  senfl,splflo,tauflo,tolnew,xrange,yc,yrange)
!
!*****************************************************************************80
!
!! TEST4 compares the sensitivities and their finite difference estimates.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparf
!
  real ( kind = 8 ) afl(ldafl,maxnfl)
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) difcof(maxcofrb)
  real ( kind = 8 ) difmax
  real ( kind = 8 ) dpar
  real ( kind = 8 ) drey
  character ( len = 2 ) eqn(neqnfl)
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gflafl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ipar
  integer ( kind = 4 ) ipivfl(maxnfl)
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jmax
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nsenfl
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) parafl(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  character ( len = 20 ) region
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) senfltmp(maxnfl,maxcofrb)
  real ( kind = 8 ) splflo(nparf)
  real ( kind = 8 ) tauflo(nparf)
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) xrange
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yrange
!
  call getsenfl(afl,area,eqn,gfl,indx,ipivfl,ldafl,maxcofrb,maxnfl,nelem, &
    neqnfl,nlband,node,np,npar,nsenfl,par,phifl,resfl,senfl)

  senfltmp(1:neqnfl,1:ncofrb) = senfl(1:neqnfl,1:ncofrb)

  dpar = drey

  call difsenfl(afl,area,difcof,dpar,eqn,gfl,gflafl,ifs, &
    ijac,indx,ipar,ipivfl,iwrite,ldafl,maxcofrb,maxelm, &
    maxnew,maxnfl, &
    ncofrb,nelem,neqnfl,nlband,node,np,npar,nparf, &
    par,parafl,phifl,region,resfl,senfl,splflo,tauflo, &
    tolnew,xrange,yc,yrange)

  imax = 0
  jmax = 0
  difmax = 0.0D+00

  do i = 1,neqnfl
    do j = 1,ncofrb
      if ( difmax <= abs(senfltmp(i,j)-senfl(i,j)) ) then
        imax = i
        jmax = j
        difmax = abs(senfltmp(i,j)-senfl(i,j))
      end if
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST4 - Results:'
  write ( *, * ) '  MAXIMUM DIFFERENCE is ',difmax
  write ( *, * ) '  I = ', imax
  write ( *, * ) '  J = ', jmax

  return
end
subroutine test5(maxcofrb,maxnfl,ncofrb,neqnfl,rb,rbase)
!
!*****************************************************************************80
!
!! TEST5 computes the product of the QR factors of the reduced basis matrix.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    30 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mhi
  integer ( kind = 4 ) mlo
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) nlo
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  real ( kind = 8 ) rbase(maxcofrb,maxcofrb)
  real ( kind = 8 ) test(10,5)
!
  ilo = 1
  ihi = min(10,neqnfl)
  jlo = 1
  jhi = min(5,ncofrb)

  do i = 1,ihi
    do j = 1,jhi
      test(i,j) = 0.0D+00
      do k = 1,ncofrb
        test(i,j) = test(i,j)+rb(i,k)*rbase(k,j)
      end do
    end do
  end do

  mhi = ihi
  mlo = ilo
  nhi = jhi
  nlo = jlo

  call prdmat(test,ihi,ilo,jhi,jlo,mhi,mlo,nhi,nlo)

  return
end
subroutine timestamp ( )
!
!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!
!  Example:
!
!    May 31 2001   9:45:54.872 AM
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
!    None
!
  implicit none
!
  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  character ( len = 8 ) date
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  character ( len = 10 )  time
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y
  character ( len = 5 ) zone
!
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
subroutine pruvpfl(gfl,indx,neqnfl,np,xc,xmax,xmin,yc,ymax,ymin)
!
!*****************************************************************************80
!
!! PRUVPFL prints the velocity and pressure.
!
!
!  Discussion:
!
!    The quantities are printed for all nodes within the user defined 
!    box (XMIN,YMIN) to (XMAX,YMAX).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 September 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) GFL(NEQNFL).
!    GFL contains the current solution estimate for the full problem,
!    containing the pressure and velocity coefficients.
!    The vector INDX must be used to index this data.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of equations (and coefficients) in the full
!    finite element system.
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    real ( kind = 8 ) XC(NP).
!    XC contains the X coordinates of the nodes.
!
!    real ( kind = 8 ) XMAX.
!    The maximum X for which a node should be displayed.
!
!    real ( kind = 8 ) XMIN.
!    The mininum X for which a node should be displayed.
! 
!    real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
!    real ( kind = 8 ) YMAX.
!    The maximum Y for which a node should be displayed.
!
!    real ( kind = 8 ) YMIN.
!    The minimum Y for which a node should be displayed.
!
  implicit none
!
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) indx(3,np)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRUVPFL - Print selected flow data'
  write ( *, * ) xmin,' = XMIN <= X <= XMAX = ',xmax
  write ( *, * ) ymin,' = YMIN <= Y <= YMAX = ',ymax
  write ( *, '(a)' ) ' '
  write ( *, * ) ' Node     X           Y      ' // &
    '      U             V             P'
  write ( *, '(a)' ) ' '

  do i = 1, np
    if ( xmin <= xc(i).and.xc(i) <= xmax.and. &
         ymin <= yc(i).and.yc(i) <= ymax) then

      i1 = indx(1,i)
      i2 = indx(2,i)
      i3 = indx(3,i)

      if ( 0 < i3 ) then
        write(*,'(i5,2g12.4,3g14.6)')i,xc(i),yc(i),gfl(i1),gfl(i2),gfl(i3)
      else
        write(*,'(i5,2g12.4,2g14.6)')i,xc(i),yc(i),gfl(i1),gfl(i2)
      end if

    end if
  end do

  return
end
subroutine pruvprb(grb,indx,maxnfl,ncofrb,nelem,node,nodelm,np,rb,xc, &
  xmax,xmin,yc,ymax,ymin)
!
!*****************************************************************************80
!
!! PRUVPRB prints the reduced velocity and pressure.
!
!
!  Discussion:
!
!    The values are printed for all nodes within the user defined box 
!    (XMIN,YMIN) to (XMAX,YMAX).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!    RB is generated by computing a finite element solution GFL,
!    which is saved for later reference as "GFLRB".
!    GFLRB is copied into the first column of RB.
!    Then, we compute the first NCOFRB derivatives of GFLRB with
!    respect to a parameter.  The first derivative
!    is stored in column 1 of RB, and so on. 
!
!    real ( kind = 8 ) XC(NP).
!    XC contains the X coordinates of the nodes.
!
!    real ( kind = 8 ) XMAX.
!    The maximum X for which a node should be displayed.
!
!    real ( kind = 8 ) XMIN.
!    The mininum X for which a node should be displayed.
!
!    real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
!    real ( kind = 8 ) YMAX.
!    The maximum Y for which a node should be displayed.
!
!    real ( kind = 8 ) YMIN.
!    The minimum Y for which a node should be displayed.
!
  implicit none
!
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nodelm(np)
  real ( kind = 8 ) prb
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  real ( kind = 8 ) urb
  real ( kind = 8 ) vrb
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xmax
  real ( kind = 8 ) xmin
  real ( kind = 8 ) xval
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) ymax
  real ( kind = 8 ) ymin
  real ( kind = 8 ) yval
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRUVPRB - Print selected flow data'
  write ( *, * ) xmin,' = XMIN <= X <= XMAX = ',xmax
  write ( *, * ) ymin,' = YMIN <= Y <= YMAX = ',ymax
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  IP   XC(IP)      YC(IP)      U(IP), V(IP), P(IP)'
  write ( *, '(a)' ) ' '

  do ip = 1, np

    if ( xmin <= xc(ip).and.xc(ip) <= xmax.and. &
         ymin <= yc(ip).and.yc(ip) <= ymax) then

      ielem = nodelm(ip)
      xval = xc(ip)
      yval = yc(ip)

      call uvprb(grb,ielem,indx,maxnfl,ncofrb,nelem,node,np,prb, &
        rb,urb,vrb,xc,xval,yc,yval)

      write(*,'(i5,2g12.4,3g14.6)')ip,xc(ip),yc(ip),urb,vrb,prb

    end if

  end do

  return
end
subroutine uvprb(grb,ielem,indx,maxnfl,ncofrb,nelem,node,np,prb, &
  rb,urb,vrb,xc,xval,yc,yval)
!
!*****************************************************************************80
!
!! UVPRB evaluates the reduced state variables at a point in an element.
!
!
!  Discusion:
!
!    The routine is given:
!
!      GRB, a set of reduced coefficients,
!      IELEM, an element,
!      XVAL, YVAL, the coordinates of a point in element IELEM,
!    and returns
!      URB, VRB, PRB, the values of the velocity and pressure defined
!      by GRB at that point.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    real ( kind = 8 ) GRB(NCOFRB).
!    GRB contains the reduced basis coefficients of the current
!    estimate of the state solution.
!
!    integer ( kind = 4 ) IELEM.
!    IELEM is the element in which the point (XVAL,YVAL) lies.
!
!    integer ( kind = 4 ) INDX(3,NP).
!    INDX(I,J) contains, for each node J, the global index of U,
!    V and P at that node, or 0 or a negative value.  The global
!    index of U, V, or P is the index of the coefficient vector
!    that contains the value of the finite element coefficient
!    associated with the corresponding basis function at the
!    given node.
!
!    integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations or coefficients allowed
!    for the full system.  MAXNFL must be used instead of NEQNFL as
!    the leading dimension of certain multi-dimensional arrays.
!
!    integer ( kind = 4 ) NCOFRB.
!    NCOFRB is the number of coefficients needed to determine
!    a particular reduced basis function.
!    NCOFRB is the sum of NBCRB and NFERB.
!
!    integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!    NELEM can be determined as 2*(NX-1)*(NY-1).
!
!    integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global index of
!    the node whose local number in J is I.
!
!    integer ( kind = 4 ) NP.
!    NP is the number of nodes used to define the finite element mesh.
!    Typically, the mesh is generated as a rectangular array, with
!    an odd number of nodes in the horizontal and vertical directions.
!    The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!    Output, real ( kind = 8 ) PRB.
!    PRB is the value of the reduced pressure at (XVAL,YVAL).
!
!    real ( kind = 8 ) RB(MAXNFL,MAXCOFRB).
!
!    RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!
!    RB is generated by computing a finite element solution GFL,
!    which is saved for later reference as "GFLRB".
!    GFLRB is copied into the first column of RB.
!    Then, we compute the first NCOFRB derivatives of GFLRB with
!    respect to a parameter.  The first derivative
!    is stored in column 1 of RB, and so on. 
!
!    Output, real ( kind = 8 ) URB.
!    URB is the value of the reduced horizontal velocity at (XVAL,YVAL).
!
!    Output, real ( kind = 8 ) VRB.
!    VRB is the value of the reduced vertical velocity at (XVAL,YVAL).
!
!    real ( kind = 8 ) XC(NP).
!    XC contains the X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XVAL.
!    XVAL is the X coordinate of the point at which the reduced
!    solution values are desired.
!
!    real ( kind = 8 ) YC(NP).
!    YC contains the Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YVAL.
!    YVAL is the Y coordinate of the point at which the reduced
!    solution values are desired.
!
  implicit none
!
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) pfl
  real ( kind = 8 ) prb
  real ( kind = 8 ) q
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  real ( kind = 8 ) ufl
  real ( kind = 8 ) urb
  real ( kind = 8 ) vfl
  real ( kind = 8 ) vrb
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xval
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yval
!
  urb = 0.0D+00
  vrb = 0.0D+00
  prb = 0.0D+00

  do i = 1,ncofrb

    ufl = 0.0D+00
    vfl = 0.0D+00
    pfl = 0.0D+00

    do j = 1,6

      ip = node(j,ielem)

      call qbf(ielem,j,w,dwdx,dwdy,nelem,node,np,xc,xval,yc,yval)

      i1 = indx(1,ip)
      ufl = ufl+rb(i1,i)*w

      i2 = indx(2,ip)
      vfl = vfl+rb(i2,i)*w

      i3 = indx(3,ip)

      if ( 0 < i3 ) then
        call bsp(q,dqdx,dqdy,ielem,j,nelem,node,np,xc,xval,yc,yval)
        pfl = pfl+rb(i3,i)*q
      end if

    end do

    urb = urb+grb(i)*ufl
    vrb = vrb+grb(i)*vfl
    prb = prb+grb(i)*pfl

  end do

  return
end
subroutine bmpcst ( costb, nparb, splbmp, taubmp, xbl, xbr, ybl, ybr )
!
!*****************************************************************************80
!
!! BMPCST evaluates the cost of the bump control.
!
!
!  Discussion:
!
!    The bump connects the points (XBL,YBL) and (XBR,YBR).
!
!    Compute its "cost" by comparing its slope to the slope of the
!    straight line that connects those two points.
!
!      COSTB = Integral (XBL <= X <= XBR) (Bump'(X) - Line'(X))**2 dX
!
!    Here, Bump(X) represents the function describing the shape
!    of the bump, and Line(X) represents the straight line which
!    simply joins the two endpoints, (XBL,YBL) and (XBR,YBR).
!
!    This integral is approximated by numerical integration.
!
!    The interval between XBL and XBR is divided into NPARB+1
!    intervals, over each of which the bump's height is described
!    by a spline.
!
!    For each such interval, pick NQUAD1 quadrature points,
!    evaluate the derivative of the bump function there, and
!    subtract the slope of the straight line.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) COSTB, the integral of the difference of the
!    derivatives of the straight line joining the two straight line
!    line segments of the bottom, and the bump that is
!    actually drawn there.  This measures the cost of bump control.
!
!    Input, integer ( kind = 4 ) NPARB, the number of parameters associated with the 
!    position and shape of the bump.
!    If NPARB = 0, the bump is replaced by a flat wall.
!
!    Input, real ( kind = 8 ) SPLBMP(NPARB+2).
!    SPLBMP contains the spline coefficients for the bump.
!
!    Input, real ( kind = 8 ) TAUBMP(NPARB+2).
!    TAUBMP contains the location of the spline abscissas for
!    the bump.  There are NPARB+2 of them, because the end values
!    of the spline are constrained to have particular values.
!
!    Input, real ( kind = 8 ) XBL, the X coordinate of the left corner
!    of the bump.
!
!    Input, real ( kind = 8 ) XBR, the X coordinate of the right corner
!    of the bump.
!
!    Input, real ( kind = 8 ) YBL, the Y coordinate of the left corner
!    of the bump.
!
!    Input, real ( kind = 8 ) YBR, the Y coordinate of the right corner
!    of the bump.
!
  implicit none
!
  integer ( kind = 4 ) nparb
!
  integer ( kind = 4 ) nquad1
  parameter (nquad1 = 5)
!
  real ( kind = 8 ) costb
  real ( kind = 8 ) cprime
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) slope
  real ( kind = 8 ) splbmp(nparb+2)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xleft
  real ( kind = 8 ) xsiquad(nquad1)
  real ( kind = 8 ) xrite
  real ( kind = 8 ) xx
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
!
  costb = 0.0D+00

  if ( nparb == 0 ) then
    return
  end if

  if ( xbr <= xbl ) then
    return
  end if
!
!  Get the Gauss weights and abscissas for one dimensional quadrature.
!
  call gquad1(nquad1,wquad1,xsiquad)
!
!  Get the slope of the line joining the endpoints of the bump.
!
  slope = (ybr-ybl) / (xbr-xbl)
!
!  Estimate the integral of the square of the difference between
!  the slope of the line and the slope of the bump over the
!  bump interval.
!
  do i = 1,nparb+1

    xleft = (dble(nparb+2-i)*xbl+dble(i-1)*xbr)/dble(nparb+1)
    xrite = (dble(nparb+1-i)*xbl+dble(i)*xbr)/dble(nparb+1)

    do j = 1,nquad1

      xx = 0.5D+00 *((1.0D+00 + xsiquad(j))*xrite+(1.0D+00 - xsiquad(j))*xleft)

      call pqdx(nparb+2,xx,taubmp,cprime,splbmp)

      costb = costb+0.5D+00 *wquad1(j)*(xrite-xleft)*(cprime-slope)**2

    end do

  end do

  return
end
subroutine bmpspl(npar,nparb,nparf,par,splbmp,taubmp,xbl,xbr,ybl,ybr)
!
!*****************************************************************************80
!
!! BMPSPL sets up or updates the spline data that describes the bump.
!
!
!  Discussion:
!
!    It does this for the target parameters and the feasible parameters.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  NPAR   Input, integer ( kind = 4 ) NPAR.
!         The number of parameters.  NPAR = NPARF + NPARB + 1.
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!         The number of parameters associated with the position and
!         shape of the bump.
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1. 
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!         PAR is the current estimate for the parameters.
!
!  SPLBMP Output, real ( kind = 8 ) SPLBMP(NPARB+2).
!         SPLBMP contains the spline coefficients for the bump.
!
!  TAUBMP Output, real ( kind = 8 ) TAUBMP(NPARB+2).
!         TAUBMP contains the location of the spline abscissas for
!         the bump.  There are NPARB+2 of them, because the end values
!         of the spline are constrained to have particular values.
!
!  XBL    Input, real ( kind = 8 ) XBL, the X coordinate of the left corner
!         of the bump.
!
!  XBR    Input, real ( kind = 8 ) XBR, the X coordinate of the right corner
!         of the bump.
!
!  YBL    Input, real ( kind = 8 ) YBL, the Y coordinate of the left corner
!         of the bump.
!
!  YBR    Input, real ( kind = 8 ) YBR, the Y coordinate of the right corner
!         of the bump.
!
  implicit none
!
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
!
  integer ( kind = 4 ) i
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) splbmp(nparb+2)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
!
  if ( nparb <= 0)return
!
!  Set up the bump arrays, including:
!
!    TAUBMP, containing the abscissas, which never change,
!    SPLBMP(I), the location of the bump at abscissa I.
!
  do i = 1,nparb+2
    taubmp(i) = ((nparb+2-i)*xbl+(i-1)*xbr)/dble(nparb+1)
  end do
!
!  Watch out!  The indexing of SPLBMP here is technically illegal.
!
    splbmp(1) = ybl
    do i = 2,nparb+1
      splbmp(i) = par(nparf+i-1)
    end do
    splbmp(nparb+2) = ybr

  return
end
subroutine bsp(q,dqdx,dqdy,ielem,iq,nelem,node,np,xc,xq,yc,yq)
!
!*****************************************************************************80
!
!! BSP evaluates the linear basis functions associated with pressure.
!
!
!  Discussion:
!
!    Here is a picture of a typical finite element associated with
!    pressure:
!
!        2
!       /|
!      / |
!     /  |
!    1---3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) Q, the value of the IQ-th basis
!    function at the point with global coordinates (XQ,YQ).
!
!  DQDX,
!  DQDY   Output, real ( kind = 8 ) DQDX, DQDY, the X and Y
!         derivatives of the IQ-th basis function at the point
!         with global coordinates (XQ,YQ).
!
!  IELEM  Input, integer ( kind = 4 ) IELEM, the global element number about which
!         we are inquiring.
!
!  IQ     Input, integer ( kind = 4 ) IQ, the index of the desired basis
!         function.  This is also the node of the reference
!         triangle which is associated with the basis function.
!
!         Basis function IQ is 1 at node IQ, and zero at the
!         other two nodes.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM).  NODE(J,I) is
!         the global node number of the J-th node in the I-th
!         element.
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  XC     Input, real ( kind = 8 ) XC(NP), the global X coordinates
!         of the element nodes.
!
!  XQ     Input, real ( kind = 8 ) XQ, the global X coordinate of
!         the point in which we are interested.
!
!    Input, real ( kind = 8 ) YC(NP), the global Y coordinates
!    of the element nodes.
!
!    Input, real ( kind = 8 ) YQ, the global Y coordinate of
!    the point in which we are interested.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) q
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) d
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iq1
  integer ( kind = 4 ) iq2
  integer ( kind = 4 ) iq3
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq
!
  if ( iq < 1 .or. 6 < iq ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'BSP - Fatal error!'
    write ( *, * ) '  The requested basis function is IQ = ',iq
    write ( *, '(a)' ) '  but only values from 1 to 6 are legal.'
    stop
  else if ( 4 <= iq .and. iq <= 6) then
    q = 0.0D+00
    dqdx = 0.0D+00
    dqdy = 0.0D+00
    return
  end if

  iq1 = iq
  iq2 = mod(iq,3)+1
  iq3 = mod(iq+1,3)+1

  i1 = node(iq1,ielem)
  i2 = node(iq2,ielem)
  i3 = node(iq3,ielem)

  d =  (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

  dqdx = (yc(i2)-yc(i3))/d
  dqdy = (xc(i3)-xc(i2))/d

  q = 1.0D+00 + dqdx*(xq-xc(i1)) + dqdy*(yq-yc(i1))

  return
end
subroutine ch_cap ( c )
!
!*****************************************************************************80
!
!! CH_CAP capitalizes a single character.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 July 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, character C, the character to capitalize.
!
  character c
  integer ( kind = 4 ) itemp
!
  itemp = ichar ( c )

  if ( 97 <= itemp .and. itemp <= 122 ) then
    c = char ( itemp - 32 )
  end if

  return
end
subroutine cavity(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb, &
  nparf,npe,nx,ny,par,region,reynld,tolnew,tolopt,tolsim,wateb,watep, &
  wateu,watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)
!
!*****************************************************************************80
!
!! CAVITY sets up the standard driven cavity problem.
!
!
!  Discussion:
!
!    This cavity has a tangential "inflow" along the top.
!
!    The strength of the top tangential flow is PAR(1).
!
!  Reference:
!
!    Janet Peterson,
!    The Reduced Basis Method for Incompressible Viscous Flow Calculations,
!    SIAM Journal of Scientific and Statistical Computing,
!    Volume 10, Number 4, pages 777-786, July 1989.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    integer ( kind = 4 ) IBS.
!    IBS is the bump shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) IBUMP.
!    IBUMP determines where isoparametric elements will be used.
!
!    0, no isoparametric elements will be used.
!       The Y coordinates of midside nodes of elements above the
!       bump will be recomputed so that the sides are straight.
!
!    1, isoparametric elements will be used only for the
!       elements which directly impinge on the bump.
!       Midside nodes of nonisoparametric elements above the
!       bump will be recomputed so that the sides are straight.
!
!    2, isoparametric elements will be used for all elements
!       which are above the bump.  All nodes above the bump
!       will be equally spaced in the Y direction.
!
!    3, isoparametric elements will be used for all elements.
!       All nodes above the bump will be equally spaced in
!       the Y direction.
!
!    integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) IOPT(MAXPAR).
!    IOPT is used during an optimization.  For each parameter I,
!    the meaning of IOPT(I) is:
!    0, the parameter value must remain fixed;
!    1, the parameter value may be varied.
!
!    integer ( kind = 4 ) MAXOPT.
!    MAXOPT is the maximum number of optimization steps.
!
!    integer ( kind = 4 ) MAXPAR.
!    MAXPAR is the maximum number of parameters allowed.
!    MAXPAR = MAXPARF + MAXPARB + 1.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
! 
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
! 
!    integer ( kind = 4 ) NPARB.
!    NPARB is the number of parameters associated with the position and
!    shape of the bump.
!
!    Note that if NPARB = 0, the bump is replaced by a flat wall.
! 
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    integer ( kind = 4 ) NPE.
!    NPE is the number of nodes per element.
! 
!    integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!
!    The number of elements along a line in the X direction is
!    NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!    integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!
!    The number of elements along a line in the Y direction is
!    NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the current value of the Reynolds number.
!    Normally, REYNLD is stored as PARA(NPARF+NPARB+1).
!
!    real ( kind = 8 ) TOLNEW.
!    TOLNEW is the convergence tolerance for the Newton iteration.
!
!    real ( kind = 8 ) TOLOPT.
!    TOLOPT is the convergence tolerance for the optimization.
!
!    real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the Picard iteration.
!
!    real ( kind = 8 ) WATEB.
!    WATEB is the multiplier of the bump control cost used
!    when computing the total cost.
!
!    real ( kind = 8 ) WATEP, WATEU, WATEV.
!
!    WATEP, WATEU and WATEV are weights used in computing the
!    cost function based on the costs of the flow discrepancy.
!
!    real ( kind = 8 ) XBL.
!    XBL is the X coordinate of the left corner of the bump.
!  
!    real ( kind = 8 ) XBR.
!    XBR is the X coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) XPROF.
!    XPROF is the X coordinate at which the profile is measured. 
!    XPROF should be a grid value!
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YBL.
!    YBL is the Y coordinate of the left corner of the bump.
!
!    real ( kind = 8 ) YBR.
!    YBR is the Y coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxpar
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(maxpar)
  character ( len = 20 ) region
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yrange
!
  ibs = 0
  ibump = 0
!
!  The "inflow" is modeled by a piecewise constant function.
!
  ifs = 0
  maxopt = 15
  nbcrb = 1
  nparb = 0
!
!  For our piecewise constant function, we specify one value.
!
  nparf = 1
  npe = 6
!
!  Peterson used a nonuniform mesh with NX = NY=25.
!
  nx = 11
  ny = 11
  region = 'cavity'
  tolnew = 0.0000000001D+00
  tolopt = 0.000000001D+00
  tolsim = 0.0000000001D+00
  wateb = 0.0D+00
  wateu = 1.0D+00
  watev = 1.0D+00
  watep = 0.0D+00
  xbl = 0.0D+00
  xbr = 0.0D+00
  xprof = 0.50
  xrange = 1.0D+00
  ybl = 0.0D+00
  ybr = 0.0D+00
  yrange = 1.0D+00
!
!  Set things that depend on other things.
!
  npar = nparf+nparb+1

  do i = 1,nparf
    iopt(i) = 0
  end do

  do i = nparf+1,nparf+nparb
    iopt(i) = 0
  end do

  iopt(nparf+nparb+1) = 1

!
!  Set the parameter that determines the tangential flow.
!
  par(1) = -1.0D+00
!
!  Set the REYNLD value.  Here, it is arbitrarily set
!  to 5.  Peterson worked with values as high as 5000.
!
  reynld = 5.0D+00
  par(2) = reynld

  return
end
subroutine cavity2(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb, &
  nparf,npe,nx,ny,par,region,reynld,tolnew,tolopt, &
  tolsim,wateb,watep,wateu,watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)
!
!*****************************************************************************80
!
!! CAVITY2 sets up the H C Lee driven cavity problem.
!
!
!  Discussion:
!
!    This cavity has a tangential "inflow" along the top, and another
!    along the bottom.
!
!    The strength of the top tangential flow is PAR(1), and the
!    strength of the bottom tangential flow is PAR(2).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    integer ( kind = 4 ) IBS.
!    IBS is the bump shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) IBUMP.
!    IBUMP determines where isoparametric elements will be used.
!
!    0, no isoparametric elements will be used.
!       The Y coordinates of midside nodes of elements above the
!       bump will be recomputed so that the sides are straight.
!
!    1, isoparametric elements will be used only for the
!       elements which directly impinge on the bump.
!       Midside nodes of nonisoparametric elements above the
!       bump will be recomputed so that the sides are straight.
!
!    2, isoparametric elements will be used for all elements
!       which are above the bump.  All nodes above the bump
!       will be equally spaced in the Y direction.
!
!    3, isoparametric elements will be used for all elements.
!       All nodes above the bump will be equally spaced in
!       the Y direction.
! 
!    integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) IOPT(MAXPAR).
!    IOPT is used during an optimization.  For each parameter I,
!    the meaning of IOPT(I) is:
!    0, the parameter value must remain fixed;
!    1, the parameter value may be varied.
!
!    integer ( kind = 4 ) MAXOPT.
!    MAXOPT is the maximum number of optimization steps.
!
!    integer ( kind = 4 ) MAXPAR.
!    MAXPAR is the maximum number of parameters allowed.
!    MAXPAR = MAXPARF + MAXPARB + 1.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NPARB.
!    NPARB is the number of parameters associated with the position and
!    shape of the bump.
!    Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    integer ( kind = 4 ) NPE.
!    NPE is the number of nodes per element.
!
!    integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!    The number of elements along a line in the X direction is
!    NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!    integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!    The number of elements along a line in the Y direction is
!    NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottom with tangential velocity specifications there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the current value of the Reynolds number.
!    Normally, REYNLD is stored as PARA(NPARF+NPARB+1).
!
!    real ( kind = 8 ) TOLNEW.
!    TOLNEW is the convergence tolerance for the Newton iteration.
!
!    real ( kind = 8 ) TOLOPT.
!    TOLOPT is the convergence tolerance for the optimization.
!
!    real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the Picard iteration.
!
!    real ( kind = 8 ) WATEB.
!    WATEB is the multiplier of the bump control cost used
!    when computing the total cost.
!
!    real ( kind = 8 ) WATEP, WATEU, WATEV.
!    WATEP, WATEU and WATEV are weights used in computing the
!    cost function based on the costs of the flow discrepancy.
!
!    real ( kind = 8 ) XBL.
!    XBL is the X coordinate of the left corner of the bump.
!
!    real ( kind = 8 ) XBR.
!    XBR is the X coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) XPROF.
!    XPROF is the X coordinate at which the profile is measured. 
!    XPROF should be a grid value!
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!    real ( kind = 8 ) YBL.
!    YBL is the Y coordinate of the left corner of the bump.
!
!    real ( kind = 8 ) YBR.
!    YBR is the Y coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxpar
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(maxpar)
  character ( len = 20 ) region
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yrange
!
  ibs = 0
  ibump = 0
!
!  The "inflow" is modeled by a piecewise constant function.
!
  ifs = 0
  maxopt = 15
  nbcrb = 1
  nparb = 0
!
!  For our piecewise constant boundary flow functions, we specify one value each,
!  top and bottom.
!
  nparf = 2
  npe = 6
!
!  Peterson used a nonuniform mesh with NX = NY=25.
!
  nx = 11
  ny = 11
  region = 'cavity2'
  tolnew = 0.0000000001
  tolopt = 0.000000001
  tolsim = 0.0000000001
  wateb = 0.0D+00
  wateu = 1.0D+00
  watev = 1.0D+00
  watep = 0.0D+00
  xbl = 0.0D+00
  xbr = 0.0D+00
  xprof = 0.50
  xrange = 1.0D+00
  ybl = 0.0D+00
  ybr = 0.0D+00
  yrange = 1.0D+00
!
!  Set things that depend on other things.
!
  npar = nparf+nparb+1

  do i = 1,nparf
    iopt(i) = 0
  end do

  do i = nparf+1,nparf+nparb
    iopt(i) = 0
  end do

  iopt(nparf+nparb+1) = 1
!
!  Set the parameters that determine the tangential flows.
!
  par(1) = -1.0D+00
  par(2) = -1.0D+00
!
!  Set the REYNLD value.  Here, it is arbitrarily set
!  to 5.  Peterson worked with values as high as 5000.
!
  reynld = 5.0D+00
  par(3) = reynld

  return
end
subroutine channl(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb, &
  nparf,npe,nx,ny,par,region,reynld,tolnew,tolopt,tolsim,wateb,watep, &
  wateu,watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)
!
!*****************************************************************************80
!
!! CHANNL sets up the standard channel problem.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    integer ( kind = 4 ) IBS.
!    IBS is the bump shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) IBUMP.
!    IBUMP determines where isoparametric elements will be used.
!    0, no isoparametric elements will be used.
!       The Y coordinates of midside nodes of elements above the
!       bump will be recomputed so that the sides are straight.
!    1, isoparametric elements will be used only for the
!       elements which directly impinge on the bump.
!       Midside nodes of nonisoparametric elements above the
!       bump will be recomputed so that the sides are straight.
!    2, isoparametric elements will be used for all elements
!       which are above the bump.  All nodes above the bump
!       will be equally spaced in the Y direction.
!    3, isoparametric elements will be used for all elements.
!       All nodes above the bump will be equally spaced in
!       the Y direction.
! 
!    integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!    integer ( kind = 4 ) IOPT(MAXPAR).
!    IOPT is used during an optimization.  For each parameter I,
!    the meaning of IOPT(I) is:
!    0, the parameter value must remain fixed;
!    1, the parameter value may be varied.
!
!    integer ( kind = 4 ) MAXOPT.
!    MAXOPT is the maximum number of optimization steps.
!
!    integer ( kind = 4 ) MAXPAR.
!    MAXPAR is the maximum number of parameters allowed.
!    MAXPAR = MAXPARF + MAXPARB + 1.
!
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
! 
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!    integer ( kind = 4 ) NPARB.
!    NPARB is the number of parameters associated with the position and
!    shape of the bump.
!
!    Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!    integer ( kind = 4 ) NPE.
!    NPE is the number of nodes per element.
! 
!    integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!    The number of elements along a line in the X direction is
!    NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!    integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!    The number of elements along a line in the Y direction is
!    NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!    character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the current value of the Reynolds number.
!    Normally, REYNLD is stored as PARA(NPARF+NPARB+1).
!
!    real ( kind = 8 ) TOLNEW.
!    TOLNEW is the convergence tolerance for the Newton iteration.
!
!    real ( kind = 8 ) TOLOPT.
!    TOLOPT is the convergence tolerance for the optimization.
!
!    real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the Picard iteration.
!
!    real ( kind = 8 ) WATEB.
!    WATEB is the multiplier of the bump control cost used
!    when computing the total cost.
!
!    real ( kind = 8 ) WATEP, WATEU, WATEV.
!    WATEP, WATEU and WATEV are weights used in computing the
!    cost function based on the costs of the flow discrepancy.
!
!    real ( kind = 8 ) XBL.
!    XBL is the X coordinate of the left corner of the bump.
!
!    real ( kind = 8 ) XBR.
!    XBR is the X coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) XPROF.
!    XPROF is the X coordinate at which the profile is measured. 
!    XPROF should be a grid value!
!
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
! 
!    real ( kind = 8 ) YBL.
!    YBL is the Y coordinate of the left corner of the bump.
!
!    real ( kind = 8 ) YBR.
!    YBR is the Y coordinate of the right corner of the bump.
!
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxpar
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(maxpar)
  character ( len = 20 ) region
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yrange
!
  ibs = 2
  ibump = 2
  ifs = 2
  maxopt = 10
  nbcrb = 1
  nparb = 3
  nparf = 1
  npe = 6
  nx = 11
  ny = 4
  region = 'channel'
  tolnew = 0.0000000001
  tolopt = 0.000000001
  tolsim = 0.0000000001
  wateb = 0.0D+00
  wateu = 1.0D+00
  watev = 1.0D+00
  watep = 0.0D+00
  xbl = 1.0D+00
  xbr = 3.0D+00
  xprof = 3.0D+00
  xrange = 10.0D+00
  ybl = 0.0D+00
  ybr = 0.0D+00
  yrange = 3.0D+00
!
!  Set things that depend on other things.
!
  npar = nparf+nparb+1

  do i = 1,nparf
    iopt(i) = 1
  end do

  do i = nparf+1,nparf+nparb
    iopt(i) = 1
  end do

  iopt(nparf+nparb+1) = 1

  par(1) = 0.5
  par(2) = 0.375
  par(3) = 0.5
  par(4) = 0.375
  reynld = 1.0D+00
  par(5) = reynld

  return
end
subroutine chrctd(string,dval,ierror,lchar)
!
!*****************************************************************************80
!
!! CHRCTD accepts a string of characters, and tries to extract a
!  real ( kind = 8 ) real number from the initial part of the
!  string.
!
!  CHRCTD will read as many characters as possible until it reaches
!  the end of the string, or encounters a character which cannot be
!  part of the number.
!
!  Legal input is:
!
!     1 blanks,
!     2 '+' or '-' sign,
!     3 integer part,
!     4 decimal point,
!     5 fraction part,
!     6 'E' or 'e' or 'D' or 'd', exponent marker,
!     7 exponent sign,
!     8 exponent integer part,
!     9 exponent decimal point,
!    10 exponent fraction part,
!    11 blanks,
!    12 final comma,
!
!  with most quantities optional.
!
!  Example:
!
!    STRING            DVAL
!
!    '1'               1.0D+00
!    '     1   '       1.0D+00
!    '1A'              1.0D+00
!    '12,34,56'        12.0D+00
!    '  34 7'          34.0D+00
!    '-1E2ABCD'        -100.0D+00
!    '-1X2ABCD'        -1.0D+00
!    ' 2E-1'           0.2
!    '23.45'           23.45
!    '-4.2E+2'         -420.0D+00
!    '17d2'            1700.0D+00
!    '-14e-2'         -0.14
!    'e2'              100.0D+00
!    '-12.73e-9.23'   -12.73 * 10.0D+00 **(-9.23)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRING Input, character ( len = * ) STRING, the string containing the
!         data to be read.  Reading will begin at position 1 and
!         terminate at the end of the string, or when no more
!         characters can be read to form a legal real.  Blanks,
!         commas, or other nonnumeric data will, in particular,
!         cause the conversion to halt.
!
!  DVAL   Output, real ( kind = 8 ) DVAL, the value that was read
!         from the string.
!
!  IERROR Output, integer ( kind = 4 ) IERROR, error flag.
!
!         0, no errors occurred.
!
!         1, 2, 6 or 7, the input number was garbled.  The
!         value of IERROR is the last type of input successfully
!         read.  For instance, 1 means initial blanks, 2 means
!         a plus or minus sign, and so on.
!
!  LCHAR  Output, integer ( kind = 4 ) LCHAR, the number of characters read from
!         STRING to form the number, including any terminating
!         characters such as a trailing comma or blanks.
!
  implicit none
!
  character chrtmp
  real ( kind = 8 ) dval
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) jbot
  integer ( kind = 4 ) jsgn
  integer ( kind = 4 ) jtop
  integer ( kind = 4 ) lchar
  logical s_eqi
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  real ( kind = 8 ) rbot
  real ( kind = 8 ) rexp
  real ( kind = 8 ) rtop
  character ( len = * ) string
!
  nchar = len(string)

  ierror = 0
  dval = 0.0D+00
  lchar = -1
  isgn = 1
  rtop = 0.0D+00
  rbot = 1.0D+00
  jsgn = 1
  jtop = 0
  jbot = 1
  ihave = 1
  iterm = 0

10    continue
  lchar = lchar+1
  chrtmp = string(lchar+1:lchar+1)
!
!  Blank character.
!
  if ( chrtmp == ' ') then
    if ( ihave == 2 .or. ihave .eq. 6 .or. ihave .eq. 7 ) then
      iterm = 1
    else if ( 1 < ihave ) then
      ihave = 11
    end if
!
!  Comma
!
  else if ( chrtmp == ',') then
    if ( ihave /= 1) then
      iterm = 1
      ihave = 12
      lchar = lchar+1
    end if
!
!  Minus sign.
!
  else if ( chrtmp == '-') then
    if ( ihave == 1) then
      ihave = 2
      isgn = -1
    else if ( ihave == 6) then
      ihave = 7
      jsgn = -1
    else
      iterm = 1
    end if
!
!  Plus sign.
!
  else if ( chrtmp == '+') then
    if ( ihave == 1) then
      ihave = 2
    else if ( ihave == 6) then
      ihave = 7
    else
      iterm = 1
    end if
!
!  Decimal point.
!
  else if ( chrtmp == '.') then
    if ( ihave < 4) then
      ihave = 4
    else if ( 6 <= ihave .and. ihave <= 8 ) then
      ihave = 9
    else
      iterm = 1
    end if
!
!  Exponent marker.
!
  else if ( s_eqi ( chrtmp,'e').or.s_eqi ( chrtmp,'d') ) then
    if ( ihave < 6) then
      ihave = 6
    else
      iterm = 1
    end if
!
!  Digit.
!
  else if ( ihave < 11.and.lge(chrtmp,'0').and.lle(chrtmp,'9') ) then

    if ( ihave <= 2) then
      ihave = 3
    else if ( ihave == 4) then
      ihave = 5
    else if ( ihave == 6.or.ihave.eq.7) then
      ihave = 8
    else if ( ihave == 9) then
      ihave = 10
    end if

    read(chrtmp,'(i1)')ndig

    if ( ihave == 3) then
      rtop = 10*rtop+ndig
    else if ( ihave == 5) then
      rtop = 10*rtop+ndig
      rbot = 10*rbot
    else if ( ihave == 8) then
      jtop = 10*jtop+ndig
    else if ( ihave == 10) then
      jtop = 10*jtop+ndig
      jbot = 10*jbot
    end if
!
!  Anything else is regarded as a terminator.
!
  else
    iterm = 1
  end if
!
!  If we haven't seen a terminator, and we haven't examined the
!  entire string, go get the next character.
!
  if ( iterm /= 1.and.lchar+1 < nchar)go to 10
!
!  If we haven't seen a terminator, and we have examined the
!  entire string, then we're done, and LCHAR is equal to NCHAR.
!
  if ( iterm /= 1.and.lchar+1 == nchar)lchar=nchar
!
!  Number seems to have terminated.  Have we got a legal number?
!  Not if we terminated in states 1, 2, 6 or 7!
!
  if ( ihave == 1.or.ihave.eq.2.or.ihave.eq.6.or.ihave.eq.7) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHRCTD - Fatal error!'
    write ( *, '(a)' ) '  Illegal or nonnumeric input!'
    return
  end if
!
!  Number seems OK.  Form it.
!
  if ( jtop == 0) then
    rexp = 1.0D+00
  else
    if ( jbot == 1) then
      rexp = 10.0D+00 **(jsgn*jtop)
    else
      rexp = dble(jsgn*jtop)
      rexp = rexp/dble(jbot)
      rexp = 10.0D+00 **rexp
    end if
  end if

  dval = dble(isgn)*rexp*rtop/rbot

  return
end
subroutine chrcti(string,intval,ierror,lchar)
!
!*****************************************************************************80
!
!! CHRCTI accepts a STRING of characters and reads an integer
!  from STRING into INTVAL.  The STRING must begin with an integer
!  but that may be followed by other information.
!
!  CHRCTI will read as many characters as possible until it reaches
!  the end of the STRING, or encounters a character which cannot be
!  part of the number.
!
!  Legal input is
!
!    blanks,
!    initial sign,
!    integer ( kind = 4 ) part,
!    blanks,
!    final comma,
!
!  with most quantities optional.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRING Input, character ( len = * ) STRING, the string containing the
!         data to be read.  Reading will begin at position 1 and
!         terminate at the end of the string, or when no more
!         characters can be read to form a legal integer.  Blanks,
!         commas, or other nonnumeric data will, in particular,
!         cause the conversion to halt.
!
!         Sample results:
!
!         STRING            INTVAL
!
!         '1'               1
!         '     1   '       1
!         '1A'              1
!         '12,34,56'        12
!         '  34 7'          34
!         '-1E2ABCD'        -100
!         '-1X2ABCD'        -1
!         ' 2E-1'           0
!         '23.45'           23
!
!  INTVAL Output, integer ( kind = 4 ) INTVAL, the integer read from the string.
!
!  IERROR Output, integer ( kind = 4 ) IERROR, error flag.
!         0 if no errors,
!         Value of IHAVE when error occurred otherwise.
!
!  LCHAR  Output, integer ( kind = 4 ) LCHAR, number of characters read from
!         STRING to form the number.
!
  implicit none
!
  character chrtmp
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ihave
  integer ( kind = 4 ) intval
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) iterm
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) lchar
  integer ( kind = 4 ) nchar
  integer ( kind = 4 ) ndig
  character ( len = * ) string
!
  nchar = len(string)

  ierror = 0
  intval = 0
  lchar = -1
  isgn = 1
  itop = 0
  ihave = 1
  iterm = 0
10    continue
  lchar = lchar+1
  chrtmp = string(lchar+1:lchar+1)

  if ( chrtmp == ' ') then
    if ( ihave == 2) then
      iterm = 1
    else if ( ihave == 3) then
      ihave = 11
    end if
  else if ( chrtmp == ',') then
    if ( ihave /= 1) then
      iterm = 1
      ihave = 12
      lchar = lchar+1
    end if
  else if ( chrtmp == '-') then
    if ( ihave == 1) then
      ihave = 2
      isgn = -1
    else
      iterm = 1
    end if
  else if ( chrtmp == '+') then
    if ( ihave == 1) then
      ihave = 2
    else
      iterm = 1
    end if
  else if ( lge(chrtmp,'0').and.lle(chrtmp,'9').and.ihave < 11) then
    ihave = 3
    read(chrtmp,'(i1)')ndig
    itop = 10*itop+ndig
  else
    iterm = 1
  end if

  if ( iterm /= 1.and.lchar+1 < nchar)go to 10
  if ( iterm /= 1.and.lchar+1 == nchar)lchar=nchar
!
!  Number seems to have terminated.  Have we got a legal number?
!
  if ( ihave == 1.or.ihave.eq.2) then
    ierror = ihave
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'CHRCTI - Fatal error!'
    write ( *, * ) '  IERROR = ',ierror
    write ( *, '(a)' ) '  Illegal or nonnumeric input:'
    write(*,'(a)')string
    return
  end if
!
!  Number seems OK.  Form it.
!
  intval = isgn*itop
  return
end
subroutine chrdb1(string)
!
!*****************************************************************************80
!
!! CHRDB1 accepts a string of characters and removes all
!  blanks and nulls, left justifying the remainder and padding with
!  blanks.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRING Input/output, character ( len = * ) STRING, the string to be
!         transformed.
!
  implicit none
!
  character chrtmp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) nchar
  character ( len = * ) string
!
  nchar = len(string)

  j = 0

  do i = 1,nchar

    chrtmp = string(i:i)
    string(i:i) = ' '

    if ( chrtmp /= ' '.and.chrtmp.ne.char(0)) then
      j = j+1
      string(j:j) = chrtmp
    end if

  end do

  return
end
subroutine chrup2(string,strng2,strng3)
!
!*****************************************************************************80
!
!! CHRUP2 copies STRING into STRNG2, up to, but not including, the
!  first occurrence of the string STRNG3.  Setting STRING = 'ABCDEFGH'
!  and STRNG3 = 'EF' results in STRNG2='ABCD'.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRING Input, character ( len = * ) STRING, the string to be copied.
!
!  STRNG2 Output, character ( len = * ) STRNG2, the copied portion of
!         STRING.
!
!  STRNG3 Input, character ( len = * ) STRNG3, the 'flag' string at which
!         the copy stops.
!
  implicit none
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) len3
  character ( len = * ) string
  character ( len = * ) strng2
  character ( len = * ) strng3
!
  len1 = len(string)
  len2 = len(strng2)
  len3 = len(strng3)

  strng2 = ' '
  i = 0
10    continue
  i = i+1
  if ( len1 < i ) then
    return
  end if

  if ( len2 < i ) then
    return
  end if

  if ( i+len3-1 <= len1) then
    if ( string(i:i+len3-1) == strng3)return
  end if

  strng2(i:i) = string(i:i)
  go to 10

end
subroutine ddetfl(afl,detlog,detman,ipivfl,lda,neqnfl,ml,mu)
!
!*****************************************************************************80
!
!! DDETFL computes the determinant of a matrix factored by DFACFL.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AFL(LDA,N).
!    AFL contains the matrix as factored by DFACFL.
!
!    Output, real ( kind = 8 ) DETLOG.
!    DETLOG is the integer part of the log base 10 of the determinant
!    of the matrix.
!
!    Output, real ( kind = 8 ) DETMAN.
!    DETMAN is the mantissa of the determinant of the matrix.
!    det(AFL) = DETMAN * 10 ** DETLOG.
!
!    Output, integer ( kind = 4 ) IPIVFL(NEQNFL), the pivot vector.
!
!    Input, integer ( kind = 4 ) LDA.
!    LDA is the leading dimension of AFL.
!    LDA must be at least 2*ML+MU+1.
!
!    Input, integer ( kind = 4 ) NEQNFL, the order of the original matrix.
!
!    Input, integer ( kind = 4 ) ML.
!    The number of diagonals below the main diagonal.
!    ML must be at least 0, and no greater than NEQNFL.
!
!    Input, integer ( kind = 4 ) MU.
!    The number of diagonals above the main diagonal.
!    MU must be at least 0, and no greater than NEQNFL.
!
  implicit none
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) neqnfl
!
  real ( kind = 8 ) afl(lda,neqnfl)
  real ( kind = 8 ) detlog
  real ( kind = 8 ) detman
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipivfl(neqnfl)
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
!
  detlog = 0.0D+00
  detman = 1.0D+00
!
  do i = 1,neqnfl

    detman = detman*afl(ml+mu+1,i)

10      continue

    if ( 10.0E+00 < abs ( detman ) ) then
      detman = detman/10.0D+00
      detlog = detlog+1.0D+00
      go to 10
    end if

20      continue

    if ( detman /= 0.0.and.abs(detman) < 1.0D+00 ) then
      detman = detman*10.0D+00
      detlog = detlog-1.0D+00
      go to 20
    end if

  end do

  do i = 1,neqnfl
    if ( ipivfl(i) /= i) then
      detman = -detman
    end if
  end do

  return
end
subroutine ddetrb(arb,detlog,detman,ipivrb,maxcofrb,ncofrb)
!
!*****************************************************************************80
!
!! DDETRB computes the determinant of a matrix factored by DFACRB.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    13 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  ARB    Input, real ( kind = 8 ) ARB(MAXNRB,NCOFRB).
!         ARB contains the matrix as factored by DFACFL.
!
!  DETLOG Output, real ( kind = 8 ) DETLOG.
!         DETLOG is the integer part of the log base 10 of the determinant
!         of the matrix.
!
!  DETMAN Output, real ( kind = 8 ) DETMAN
!         DETMAN is the mantissa of the determinant of the matrix.
!         det(ARB) = DETMAN * 10 ** DETLOG.
!
!  IPIVRB Input, integer ( kind = 4 ) IPIVRB(NCOFRB).
!         The pivot vector.
!
!  MAXNRB Input, integer ( kind = 4 ) MAXNRB.
!         The leading dimension of the array ARB.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB.
!         The order of the original matrix.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) ncofrb
!
  real ( kind = 8 ) arb(maxcofrb,ncofrb)
  real ( kind = 8 ) detlog
  real ( kind = 8 ) detman
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipivrb(ncofrb)
!
  detlog = 0.0D+00
  detman = 1.0D+00
!
  do i = 1,ncofrb

    detman = detman*arb(i,i)

10      continue

    if ( 10.0E+00 < abs(detman) ) then
      detman = detman/10.0D+00
      detlog = detlog+1.0D+00
      go to 10
    end if

20      continue

    if ( detman /= 0.0.and.abs(detman) < 1.0D+00 ) then
      detman = detman*10.0D+00
      detlog = detlog-1.0D+00
      go to 20
    end if

  end do

  do i = 1,ncofrb
    if ( ipivrb(i) /= i) then
      detman = -detman
    end if
  end do

  return
end
subroutine delhms ( time1, time2, nsec )
!
!*****************************************************************************80
!
!! DELHMS returns the number of seconds between TIME1 and TIME2.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    04 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character*10 TIME1, TIME2, two times, in decimal seconds,
!    represented as character strings.
!
!    Output, integer ( kind = 4 ) NSEC, the number of elapsed seconds.
!
  integer ( kind = 4 ) nsec
  real rsec1
  real rsec2
  character ( len = 10 ) time1
  character ( len = 10 ) time2
!
  read ( time1, '(f10.3)' ) rsec1
  read ( time2, '(f10.3)' ) rsec2

  nsec = int ( rsec2 - rsec2 )

  return
end
subroutine dfacfl(afl,lda,n,ml,mu,ipivfl,info)
!
!*****************************************************************************80
!
!! DFACFL factors a real ( kind = 8 ) band matrix by elimination.
!
!
!  Discussion:
!
!    DFACFL is a simplified version of the LINPACK routine DGBFA.
!
!    In order to use DFACFL, it is necessary to store the matrix AFL
!    in "LINPACK General Band Storage" format.
!
!    If AFL is a band matrix, the following program segment
!    will set up the compressed matrix properly:
!
!        m = ml+mu+1
!        do j = 1,n
!          i1 = max(1,j-mu)
!          i2 = min(n,j+ml)
!          do i = i1,i2
!            k = i-j+m
!            afl(k,j) = Entry I, J
!          end do
!        end do
!
!    This uses rows ML+1 through 2*ML+MU+1 of the array AFL.
!    In addition, the first ML rows in ABD are used for
!    elements generated during the triangularization because of pivoting.
!    The total number of rows needed in AFL is 2*ML+MU+1.
!    The ML+MU by ML+MU upper left triangle and the
!    ML by ML lower right triangle are not referenced.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    22 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  AFL    Input/output, real ( kind = 8 ) AFL(LDA,N).
!         On input, AFL contains the matrix in band storage.  The
!         columns of the matrix are stored in the columns of AFL and
!         the diagonals of the matrix are stored in rows
!         ML+1 through 2*ML+MU+1 of AFL.
!
!         On output, an upper triangular matrix in band storage and
!         the multipliers which were used to obtain it.
!         The factorization can be written AFL = L*U where
!         L is a product of permutation and unit lower
!         triangular matrices and U is upper triangular.
!
!  LDA    Input, integer ( kind = 4 ) LDA.
!         The leading dimension of the array AFL.
!         LDA must be at least 2*ML+MU+1.
!
!  N      Input, integer ( kind = 4 ) N.
!         The order of the original matrix.
!
!  ML     Input, integer ( kind = 4 ) ML.
!         The number of diagonals below the main diagonal.
!         ML must be at least 0, and no greater than N.
!
!  MU     Input, integer ( kind = 4 ) MU.
!         The number of diagonals above the main diagonal.
!         MU must be at least 0, and no greater than N.
!
!  IPIVFL Output, integer ( kind = 4 ) IPIVFL(N).
!         An integer vector of pivot indices needed by DSOLFL.
!
!  INFO   Output, integer ( kind = 4 ) INFO.
!         = 0  normal value.
!         = K  if U(K,K)  ==  0.0.  In this case, the matrix is exactly
!         numerically singular, and DSOLFL should not be called to attempt
!         a linear solution.
!
  implicit none
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) afl(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i0
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivfl(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) ju
  integer ( kind = 4 ) jz
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mm
  integer ( kind = 4 ) mu
  real ( kind = 8 ) t
!
  m = ml+mu+1
  info = 0
!
!  Zero out the initial fill-in columns of the matrix.
!
  j1 = min(n,m)-1
  do jz = mu+2,j1
    i0 = m+1-jz
    do i = i0,ml
      afl(i,jz) = 0.0D+00
    end do
  end do

  jz = j1
  ju = 0
!
!  Carry out Gaussian elimination with partial pivoting
!
  do k = 1,n-1
!
!  Zero out the next fill-in column.
!
    jz = jz+1
    if ( jz <= n) then
      do i = 1,ml
         afl(i,jz) = 0.0D+00
      end do
    end if
!
!  Find L = pivot index
!
    lm = min(ml,n-k)

    l = m
    do i = m+1,m+lm
      if ( abs(afl(l,k)) < abs(afl(i,k)) ) then
        l = i 
      end if
    end do

    ipivfl(k) = l+k-m
!
!  A zero pivot means the matrix is singular.
!
    if ( afl(l,k) == 0.0D+00 ) then
      info = k
    else
!
!  Interchange rows unless the pivot row is already on the diagonal.
!
      if ( l /= m) then
        t = afl(l,k)
        afl(l,k) = afl(m,k)
        afl(m,k) = t
      end if
!
!  Compute the multipliers that form the lower diagonal entries of
!  the L factor.
!
      do i = m+1,m+lm
        afl(i,k) = -afl(i,k)/afl(m,k)
      end do
!
!  Row elimination with column indexing.
!
      ju = max(ju,mu+ipivfl(k))
      ju = min(ju,n)
      mm = m

      do j = k+1,ju
        l = l-1
        mm = mm-1
       
        t = afl(l,j)
        if ( l /= mm) then
          afl(l,j) = afl(mm,j)
          afl(mm,j) = t
        end if
       
        do i = 1,lm
          afl(mm+i,j) = afl(mm+i,j)+afl(m+i,k)*t
        end do

      end do

    end if
   
  end do
 
  ipivfl(n) = n
 
  if ( afl(m,n) == 0.0D+00 ) then
    info = n
  end if
 
  return
end
subroutine dfacrb(a,lda,n,ipivot,info)
!
!*****************************************************************************80
!
!! DFACRB factors a real ( kind = 8 ) dense matrix.
!
!
!  Discussion:
!
!    DFACRB is similar to the LINPACK routine DGEFA, but does not call
!    any subroutines or functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1996.
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!
!    On input, A contains the N by N matrix to be factored.
!
!    On output, A contains the L and U factors of the matrix, in
!    compressed storage.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!    LDA may be larger than N, but must not be smaller than N.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Output, integer ( kind = 4 ) IPIVOT(N), the pivot array.
!
!    Output, integer ( kind = 4 ) INFO, an error flag.
!
!    INFO = 0, no error, the matrix was factored.
!    INFO = K, the K-th pivot U(K,K) was zero.
!
  implicit none
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
  info = 0

  do k = 1,n-1
!
!  Find the pivot row L.
!
    l = k
    do i = k+1,n
      if ( abs(a(i,l)) <  abs(a(i,k)) ) then
        l = i
      end if
    end do

    ipivot(k) = l
!
!  Check for a zero pivot.
!
    if ( a(l,k) == 0.0D+00 ) then

      info = k
      return

    end if
!
!  Check to see whether we must swap rows L and K.
!
    if ( l /= k) then
      t = a(l,k)
      a(l,k) = a(k,k)
      a(k,k) = t
    end if
!
!  Rescale the pivot row so that A(K,K) = 1.
!
    do i = k+1,n
      a(i,k) = -a(i,k)/a(k,k)
    end do
!
!  Wipe out the entries below A(K,K).
!
    do j = k+1,n

      t = a(l,j)

      if ( l /= k) then
        a(l,j) = a(k,j)
        a(k,j) = t
      end if

      do i = k+1,n
        a(i,j) = a(i,j)+t*a(i,k)
      end do

    end do

  end do

  ipivot(n) = n

  if ( a(n,n) == 0.0D+00 ) then
    info = n
  end if

  return
end
subroutine difset ( difcof, h, iwrite, ncof )
!
!*****************************************************************************80
! 
!! DIFSET computes the NCOF coefficients for a centered finite difference
!  estimate of the (NCOF-1)-th derivative of a function.
!
!
!  The estimate has the form
!
!    FDIF(NCOF-1,X) = Sum (I=1 to NCOF) COF(I) * F(X(I))
!
!  To understand the computation of the coefficients, it is enough
!  to realize that the first difference approximation is
!
!    FDIF(1,X) = F(X+DX) - F(X-DX) ) / (2*DX)
!
!  and that the second difference approximation can be regarded as
!  the first difference approximation repeated:
!
!    FDIF(2,X) = FDIF(1,X+DX) - FDIF(1,X-DX) / (2*DX)
!           = F(X+2*DX) - 2 F(X) + F(X-2*DX) / (4*DX)
!
!  and so on for higher order differences.
!
!  Thus, the next thing to consider is the integer coefficients of
!  the sampled values of F, which are clearly the Pascal coefficients,
!  but with an alternating negative sign.  In particular, if we
!  consider row I of Pascal's triangle to have entries J = 0 through I,
!  then P(I,J) = P(I-1,J-1) - P(I-1,J), where P(*,-1) is taken to be 0,
!  and P(0,0) = 1.
!
!     1
!    -1  1
!     1 -2   1
!    -1  3  -3   1
!     1 -4   6  -4   1
!    -1  5 -10  10  -5  1
!     1 -6  15 -20  15 -6 1
!
!  Next, we note that the denominator of the approximation for the
!  I-th derivative will be (2*DX)**I.
!
!  And finally, we must consider the location of the NDIF sampling
!  points for F:
!
!    X-NDIF*DX, X-(NDIF-2)*DX, X-(NDIF-4)*DX, ...,
!    X+(NDIF-4)*DX, X+(NDIF-2*DX), X+(NDIF-1)*DX.
!
!
!  Thus, a formula for evaluating FDIF(NDIF,X) is
!
!    fdif = 0.0D+00
!    ncof = ndif+1
!    do i = 1,ncof
!      xi = x+(2*(i-1)-ncof)*h
!      fdif = fdif+difcof(i)*f(xi)
!    end do
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DIFCOF Output, real DIFCOF(NCOF), the coefficients needed to approximate
!         the (NCOF-1)-th derivative of a function F.
!
!  H      Input, real H, the half spacing between points.  H must
!         be positive.
!
!  NCOF   Input, integer ( kind = 4 ) NCOF.
!         NCOF is the number of coefficients desired, which also
!         determines NDIF = NCOF-1, the derivative being estimated.
!
  implicit none
!
  integer ( kind = 4 ) ncof
!
  real ( kind = 8 ) difcof(ncof)
  real ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iwrite
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ndif
!
  ndif = ncof-1

  if ( ndif < 0) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFSET - Fatal error!'
    write ( *, * ) '  Derivative order NDIF = ',ndif
    write ( *, '(a)' ) '  but NDIF must be at least 0.'
    stop
  end if

  if ( h <= 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFSET - Fatal error!'
    write ( *, * ) '  The half sampling spacing is H = ',H
    write ( *, '(a)' ) '  but H must be positive.'
    stop
  end if

  do i = 1,ncof

    difcof(i) = 1.0D+00

    do j = i-1,2,-1
      difcof(j) = -difcof(j)+difcof(j-1)
    end do

    if ( 1 < i ) then
      difcof(1) = -difcof(1)
    end if

  end do

  if ( 2 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFSET - Unnormalized coefficients:'
    do i = 1,ncof
      write(*,'(i6,g14.6)')i,difcof(i)
    end do
  end if

  do i = 1,ncof
    difcof(i) = difcof(i)/(2.0D+00 *h)**ndif
  end do

  if ( 2 <= iwrite ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'DIFSET - Normalized coefficients:'
    do i = 1,ncof
      write(*,'(i6,g14.6)')i,difcof(i)
    end do
  end if

  return
end
subroutine discst(costp,costu,costv,gfl,gfltar,indx,neqnfl,np,nprof,ny,yc)
!
!*****************************************************************************80
!
!! DISCST computes the discrepancy integrals for the pressure,
!  horizontal and vertical velocities, along the profile line.
!
!  Discussion:
!
!    This integration scheme assumes that the profile line, and
!    the element sides that define it, are straight.  Otherwise,
!    the integration scheme used is not correct.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    01 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  COSTP  Output, real ( kind = 8 ) COSTP.
!
!         The integral of the difference between
!         the computed and target pressure functions along the
!         profile line.
!
!  COSTU  Output, real ( kind = 8 ) COSTU.
!
!         The integral of the difference between
!         the computed and target horizontal velocity functions along
!         the profile line.
!
!  COSTV  Output, real ( kind = 8 ) COSTV.
!
!         The integral of the difference between
!         the computed and target vertical velocity functions along
!         the profile line.
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL), the current solution
!         estimate for the full problem.
!
!  GTARFL Input, real ( kind = 8 ) GTARFL(NEQNFL), the target solution vector.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of equations in the full system.
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  NPROF  Input, integer ( kind = 4 ) NPROF(2*MAXNY-1).
!
!         NPROF contains the numbers of the nodes along the profile
!         line.
!
!  NY     Input, integer ( kind = 4 ) NY.
!
!         NY controls the spacing of nodes and elements in
!         the Y direction.  There are 2*NY-1 nodes along various
!         lines in the Y direction.
!
!  YC     Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) nquad1
!
  parameter (nquad1 = 5)
!
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) ny
!
  real ( kind = 8 ) bval
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gfltar(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) npol
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) pcof(2)
  real ( kind = 8 ) pval
  real ( kind = 8 ) ucof(3)
  real ( kind = 8 ) uval
  real ( kind = 8 ) vcof(3)
  real ( kind = 8 ) vval
  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xsiquad(nquad1)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) ypol(3)
  real ( kind = 8 ) yval
!
!  Get the weights and abscissas to approximate a line integral.
!
  call gquad1(nquad1,wquad1,xsiquad)
!
!  Compute the integral of the difference squared between the
!  current velocity and the target values.
!
  costu = 0.0D+00
  costv = 0.0D+00
!
!  The line along which we integrate is broken into NY-1
!  subintervals, over each of which, U and V are represented
!  by quadratic functions.
!
  do i = 1,ny-1
!
!  Get the values of U and V at the beginning, middle, and
!  end of the subinterval.  Use these to compute the quadratic
!  representation of U and V for any point on the subinterval.
!
    ylo = yc(nprof(2*i-1))
    yhi = yc(nprof(2*i+1))

    npol = 3

    do k = 1,npol

      ii = 2*i-2+k
      ypol(k) = yc(nprof(ii))

      j = indx(1,nprof(ii))
      ucof(k) = gfl(j)-gfltar(j)

      j = indx(2,nprof(ii))
      vcof(k) = gfl(j)-gfltar(j)

    end do
!
!  Evaluate the discrepancy at each quadrature point.
!
    do j = 1,nquad1

      yval = 0.5D+00 *((1.0D+00 + xsiquad(j))*ylo+(1.0D+00 -xsiquad(j))*yhi)

      uval = 0.0D+00
      vval = 0.0D+00

      do k = 1,npol
        call lbase(k,npol,bval,ypol,yval)
        uval = uval+bval*ucof(k)
        vval = vval+bval*vcof(k)
      end do

      costu = costu+0.5D+00 *wquad1(j)*(yhi-ylo)*uval**2
      costv = costv+0.5D+00 *wquad1(j)*(yhi-ylo)*vval**2

    end do
  end do
!
!  Compute the square root of the integral of the difference
!  squared between the current pressure and the target values.
!
  costp = 0.0D+00

  do i = 1,ny-1

    ylo = yc(nprof(2*i-1))
    yhi = yc(nprof(2*i+1))

    npol = 2

    do k = 1,npol

      ii = 2*i-3+2*k

      ypol(k) = yc(nprof(ii))

      j = indx(3,nprof(ii))
      if ( j <= 0) then
        pcof(k) = 0.0D+00
      else
        pcof(k) = gfl(j)-gfltar(j)
      end if

    end do

    do j = 1,nquad1

      yval = 0.5D+00 *((1.0D+00 + xsiquad(j))*ylo+ (1.0D+00 -xsiquad(j))*yhi)

      pval = 0.0D+00

      do k = 1,npol
        call lbase(k,npol,bval,ypol,yval)
        pval = pval+bval*pcof(k)
      end do

      costp = costp+0.5D+00 *wquad1(j)*(yhi-ylo)*pval**2

    end do
  end do

  return
end
subroutine dsolfl ( afl, lda, n, ml, mu, ipivfl, b )
!
!*****************************************************************************80
!
!! DSOLFL solves a (full) banded linear system.
!
!
!  Discussion:
!
!    The linear system has the form:
!
!      AFL*X = B
!
!    where AFL, X, and B are real ( kind = 8 ), and AFL is a banded matrix
!    which has already been decomposed into LU factors by DFACFL.
!
!    DSOLFL is a simplied version of the LINPACK routine DGBSL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AFL(LDA,N).
!    The factored matrix produced by DFACFL.
!
!    Input, integer ( kind = 4 ) LDA.
!    The leading dimension of the array AFL.
!
!    Input, integer ( kind = 4 ) N.
!    The order of the original matrix.
!
!    Input, integer ( kind = 4 ) ML.
!    The number of diagonals below the main diagonal.
!
!    Input, integer ( kind = 4 ) MU.
!    The number of diagonals above the main diagonal.
!
!    Input, integer ( kind = 4 ) IPIVFL(N).
!    The pivot vector from DFACFL.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector X.
!
  implicit none
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) afl(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipivfl(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lm
  integer ( kind = 4 ) m
  integer ( kind = 4 ) ml
  integer ( kind = 4 ) mu
  real ( kind = 8 ) t
!
  m = mu+ml+1
!
!  First solve L*Y = B.
!
  if ( ml /= 0) then
 
    do k = 1,n-1
   
      lm = min(ml,n-k)
      l = ipivfl(k)

      if ( l /= k) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if
     
      do i = 1,lm
        b(k+i) = b(k+i)+afl(m+i,k)*b(k)
      end do
     
    end do
  end if
!
!  Now solve U*X = Y.
!
  do k = n,1,-1
 
    if ( afl(m,k) == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DSOLFL - Fatal error!'
      write ( *, * ) '  Pivot K = ',k,' is zero.'
      stop
    else
      b(k) = b(k)/afl(m,k)
    end if
   
    lm = min(k,m)-1
    la = m-lm
    lb = k-lm
    do i = 1,lm
      b(lb-1+i) = b(lb-1+i)-afl(la-1+i,k)*b(k)
    end do

  end do
 
  return
end
subroutine dsolrb ( a, lda, n, ipivot, b )
!
!*****************************************************************************80
!
!! DSOLRB solves a (reduced) dense linear system.
!
!
!  Discussion:
!
!    The linear system has the form:
!
!      A*X = B
!
!    where A is a full storage real ( kind = 8 ) array which has been
!    LU-factored by DFACRB.
!
!    DSOLRB is similar to the LINPACK routine DGESL, but does not call
!    any subroutines or functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(LDA,N).
!    A contains the LU factors of a matrix, as computed by DFACRB.
!
!    Input, integer ( kind = 4 ) LDA.
!    LDA is the leading dimension of the matrix A.
!
!    Input, integer ( kind = 4 ) N.
!    N is the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPIVOT(N).
!    IPIVOT is the pivot vector computed by DFACRB.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, B is the right hand side of the linear system.
!    On output, B is the solution of the linear system.
!
  implicit none
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipivot(n)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real ( kind = 8 ) t
!
!  First solve L*Y = B.
!
  do k = 1,n-1

    l = ipivot(k)

    t = b(l)
    if ( l /= k) then
      b(l) = b(k)
      b(k) = t
    end if

    do i = k+1,n
      b(i) = b(i)+t*a(i,k)
    end do

  end do
!
!  Now solve U*X = Y.
!
  do k = n,1,-1

    b(k) = b(k)/a(k,k)
    t = -b(k)

    do i = 1,k-1
      b(i) = b(i)+t*a(i,k)
    end do

  end do

  return
end
function dveq(n,dvec1,dvec2)
!
!*****************************************************************************80
!
!! DVEQ returns .TRUE. if the N elements of the real ( kind = 8 )
!  vectors DVEC1 and DVEC2 are equal.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  N      Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!  DVEC1,
!  DVEC2  Input, real ( kind = 8 ) DVEC1(N), DVEC2(N), the two vectors
!         to be compared.
!
!  DVEQ   Output, logical DVEQ.
!         DVEQ is .TRUE. if all N elements of DVEC1 and DVEC2 are equal,
!         and .FALSE. otherwise.
!
  implicit none
!
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) dvec1(n)
  real ( kind = 8 ) dvec2(n)
  logical dveq
  integer ( kind = 4 ) i
!
  dveq = .false.
 
  do i = 1,n
    if ( dvec1(i) /= dvec2(i))return
  end do
 
  dveq = .true.
 
  return
end
function dvneq(n,dvec1,dvec2)
!
!*****************************************************************************80
!
!! DVNEQ returns .TRUE. if any of the N elements of the real ( kind = 8 )
!  vectors DVEC1 and DVEC2 are not equal.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N.
!    N is the number of entries in the vectors.
!
!    Input, real ( kind = 8 ) DVEC1(N), DVEC2(N).
!    DVEC1 and DVEC2 are the two vectors to be compared.
!
!    Output, logical DVNEQ.
!    DVNEQ is .TRUE. if any elements of DVEC1 and DVEC2 differ,
!    and .FALSE. otherwise.
!
  implicit none
!
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) dvec1(n)
  real ( kind = 8 ) dvec2(n)
  logical dvneq
  integer ( kind = 4 ) i
!
  dvneq = .true.
 
  do i = 1,n
    if ( dvec1(i) /= dvec2(i))return
  end do
 
  dvneq = .false.
 
  return
end
subroutine fact(n,factn)
!
!*****************************************************************************80
!
!! FACT computes the (real) factorial of a nonnegative integer.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the nonnegative value for which N!
!    is desired.
!
!    Output, real ( kind = 8 ) FACTN, the factorial of N.
!
  implicit none
!
  real ( kind = 8 ) factn
  integer ( kind = 4 ) i
  integer ( kind = 4 ) n
!
  if ( n < 0) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'FACT - Fatal error!'
    write ( *, * ) '  Negative input argument is N = ',n
    stop
  end if

  factn = 1.0D+00
  do i = 1,n
    factn = factn*dble(i)
  end do

  return
end
subroutine getcst(cost,costb,costp,costu,costv,gfl,gfltar,indx,neqnfl,np, &
  nparb,nprof,ny,splbmp,taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
!
!*****************************************************************************80
!
!! GETCST is given the value of the solution, GFL, the target
!  solution GTARFL, and information about the shape of the bump,
!  and returns the value of the overall and individual cost
!  functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  COST   Output, real ( kind = 8 ) COST, the weighted cost.
!
!  COSTB  Output, real ( kind = 8 ) COSTB.
!
!         COSTB is the integral of the difference of the
!         derivatives of the straight line joining the two straight line
!         line segments of the bottom, and the bump that is
!         actually drawn there.
!
!         This measures the cost of bump control.
!
!  COSTP  Output, real ( kind = 8 ) COSTP.
!
!         The integral of the difference between
!         the computed and target pressure functions along the
!         profile line.
!
!  COSTU  Output, real ( kind = 8 ) COSTU.
!
!         The integral of the difference between
!         the computed and target horizontal velocity functions along
!         the profile line.
!
!  COSTV  Output, real ( kind = 8 ) COSTV.
!
!         The integral of the difference between
!         the computed and target vertical velocity functions along
!         the profile line.
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL), the current solution
!         estimate for the full problem.
!
!  GTARFL Input, real ( kind = 8 ) GTARFL(NEQNFL), the target solution vector.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of equations in the full system.
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!
!         The number of parameters associated with the position and
!         shape of the bump.
!
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPROF  Input, integer ( kind = 4 ) NPROF(2*MAXNY-1).
!
!         NPROF contains the numbers of the nodes along the profile
!         line.
!
!  NY     Input, integer ( kind = 4 ) NY.
!
!         NY controls the spacing of nodes and elements in
!         the Y direction.  There are 2*NY-1 nodes along various
!         lines in the Y direction.
!
!  SPLBMP Input, real ( kind = 8 ) SPLBMP(NPARB+2).
!
!         SPLBMP contains the spline coefficients for the bump.
!
!  TAUBMP Input, real ( kind = 8 ) TAUBMP(NPARB+2).
!
!         TAUBMP contains the location of the spline abscissas for
!         the bump.  There are NPARB+2 of them, because the end values
!         of the spline are constrained to have particular values.
!
!  WATEB  Input, real ( kind = 8 ) WATEB.
!
!         WATEB is the multiplier of the bump control cost used
!         when computing the total cost.
!
!  WATEP,
!  WATEU,
!  WATEV  Input, real ( kind = 8 ) WATEP, WATEU, WATEV.
!
!         These are weights used in computing the overall cost
!         function based on the costs of the flow discrepancy.
!
!  XBL    Input, real ( kind = 8 ) XBL, the X coordinate of the left corner
!         of the bump.
!
!  XBR    Input, real ( kind = 8 ) XBR, the X coordinate of the right corner
!         of the bump.
!
!  YBL    Input, real ( kind = 8 ) YBL, the Y coordinate of the left corner
!         of the bump.
!
!    Input, real ( kind = 8 ) YBR, the Y coordinate of the right corner
!    of the bump.
!
!    Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) ny
!
  real ( kind = 8 ) cost
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gfltar(neqnfl)
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) nprof(2*ny-1)
  real ( kind = 8 ) splbmp(nparb+2)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
!
  call bmpcst(costb,nparb,splbmp,taubmp,xbl,xbr,ybl,ybr)

  call discst(costp,costu,costv,gfl,gfltar,indx,neqnfl,np,nprof,ny,yc)

  cost = wateb*costb+watep*costp+wateu*costu+watev*costv

  return
end
subroutine gquad1(nquad1,wquad1,xsiquad)

!*****************************************************************************80
!
!! GQUAD1 defines a 1 dimensional Gauss quadrature rule.
!
!
!  Discussion:
!
!    GQUAD1 returns the weights and abscissas for a 1 dimensional,
!    3 or 5 point Gauss quadrature rule defined on the interval [-1,1].
!
!    The integral of a function F(X) over the interval [-1,1]
!
!      Integral (-1 to 1) F(X) DX
!
!    may then be approximated by
!
!      Sum (I = 1 to NQUAD1) WQUAD1(I) * F(XSIQUAD(I))
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  NQUAD1 Input, integer ( kind = 4 ) NQUAD1.
!         The user specifies the rule desired by setting NQUAD1
!         to 3 or 5.  Any other value is illegal, and will cause
!         GQUAD1 to stop.
!
!  WQUAD1 Output, real ( kind = 8 ) WQUAD1(NQUAD1).
!         WQUAD1(I) is the weight factor corresponding to the
!         I-th quadrature point.
!
!  XSIQUAD
!         Output, real ( kind = 8 ) XSIQUAD(NQUAD1).
!         XSIQUAD(I) is the I-th quadrature point.
!
  implicit none
!
  integer ( kind = 4 ) nquad1
!
  real ( kind = 8 ) wquad1(nquad1)
  real ( kind = 8 ) xsiquad(nquad1)
!
  if ( nquad1 == 3) then

    xsiquad(1) = -0.7745966692
    xsiquad(2) =  0.0D+00
    xsiquad(3) =  0.7745966692

    wquad1(1) = 5.0D+00 / 9.0D+00
    wquad1(2) = 8.0D+00 / 9.0D+00
    wquad1(3) = 5.0D+00 / 9.0D+00

  else if ( nquad1 == 5) then

    xsiquad(1) = -0.906179845938664D+00
    xsiquad(2) = -0.538469310105683D+00
    xsiquad(3) =  0.0D+00
    xsiquad(4) =  0.538469310105683D+00
    xsiquad(5) =  0.906179845938664D+00

    wquad1(1) = 0.236926885056189D+00
    wquad1(2) = 0.478628670499366D+00
    wquad1(3) = 0.568888888888889D+00
    wquad1(4) = 0.478628670499366D+00
    wquad1(5) = 0.236926885056189D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'GQuad1 - Fatal error!'
    write ( *, '(a)' ) '  An illegal value of NQUAD1 was input.'
    write ( *, '(a)' ) '  Only NQUAD1 = 3 or 5 are legal.'
    write ( *, * ) '  The input value was ',nquad1
    write ( *, '(a)' ) '  The code is stopping now.'
    stop

  end if

  return
end
subroutine grid(gridx,i,ihi,ilo,x,xhi,xlo)

!*****************************************************************************80
!
!! GRID computes the X or Y coordinate of the I-th gridpoint.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GRIDX  Input, character*20 GRIDX.
!         GRIDX tells how the finite element nodes should be layed out
!         in the X direction.
!         'uniform' makes them equally spaced.
!         'cos' uses the COS function to cluster them near edges.
!         'sqrtsin' uses the SQRT(SIN()) function to cluster near edges.
!
!  I      Input, integer ( kind = 4 ) I.
!         I is the index of the grid point whose X coordinate is to
!         be computed.  Normally, ILO <= I <= IHI.
!
!  IHI,
!  ILO    Input, integer ( kind = 4 ) IHI, ILO.
!         ILO is the index of the grid point whose X coordinate is XLO,
!         IHI is the same for XHI.
!
!  X      Output, real ( kind = 8 ) X.
!         X is the X coordinate of the I-th grid point, according to
!         the specified scheme.
!
!  XHI,
!  XLO    Input, real ( kind = 8 ) XHI, XLO.
!         XLO is the X coordinate of grid point ILO, and XHI
!         is the X coordinate of grid point IHI.
!
  implicit none
!
  real ( kind = 8 ) pi
  parameter (pi = 3.14159265)
!
  character ( len = 20 ) gridx
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  logical s_eqi
  real ( kind = 8 ) s
  real ( kind = 8 ) theta
  real ( kind = 8 ) thi
  real ( kind = 8 ) tlo
  real ( kind = 8 ) x
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
!
  if ( s_eqi ( gridx,'uniform')) then
    x = (dble(ihi-i)*xlo+dble(i-ilo)*xhi)/dble(ihi-ilo)
  else if ( s_eqi ( gridx,'sin')) then
    tlo = -pi/2.0D+00
    thi = pi/2.0D+00
    theta = (dble(ihi-i)*tlo + dble(i-ilo)*thi)/dble(ihi-ilo)
    s = sin(theta)
    x = ((1.0D+00 -s)*xlo+(s+1.0D+00 )*xhi)/2.0D+00
!
!  Equivalent to 'SIN'.
!
  else if ( s_eqi ( gridx,'cos')) then
    tlo = -pi
    thi = 0.0D+00
    theta = (dble(ihi-i)*tlo + dble(i-ilo)*thi)/dble(ihi-ilo)
    x = ((1.0D+00 -cos(theta))*xlo+(1.0D+00 + cos(theta))*xhi)/2.0D+00
  else if ( s_eqi ( gridx,'sqrtsin')) then
    tlo = -pi/2.0D+00
    thi = pi/2.0D+00
    theta = (dble(ihi-i)*tlo + dble(i-ilo)*thi)/dble(ihi-ilo)
    if ( 0.0D+00 <= sin(theta) ) then
      s = sqrt(sin(theta))
    else
      s = -sqrt(-sin(theta))
    end if
    x = ((1.0D+00 -s)*xlo+(s+1.0D+00 )*xhi)/2.0D+00
  end if

  return
end
subroutine intprs(gfl,indx,nelem,neqnfl,node,np,p)

!*****************************************************************************80
!
!! INTPRS interpolates the pressure at the midside nodes.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!         GFL is the current solution estimate for the full problem,
!         containing pressure and velocity coefficients.  The vector
!         INDX must be used to index this data.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the global index of U,
!         V and P at that node, or 0 or a negative value.  The global
!         index of U, V, or P is the index of the coefficient vector
!         that contains the value of the finite element coefficient
!         associated with the corresponding basis function at the
!         given node.
!  
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL.
!         NEQNFL is the number of equations (and coefficients) in the full
!         finite element system.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!
!         NODE(I,J) contains, for an element J, the global index of
!         the node whose local number in J is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!           Global nodes   Elements      NODE
!                                                          1  2  3  4  5  6
!           74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                          |  /   /|
!           73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                          |/   /  |
!           72  82  92     2   1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  P      Input, real P(NP), the pressure.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) in3
  integer ( kind = 4 ) in4
  integer ( kind = 4 ) in5
  integer ( kind = 4 ) in6
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p(np)
!
!  For each element,...
!
  do i = 1,nelem
!
!  Get the six global node numbers.
!
    in1 = node(1,i)
    in2 = node(2,i)
    in3 = node(3,i)
    in4 = node(4,i)
    in5 = node(5,i)
    in6 = node(6,i)
!
!  Read off the three computed values, and average the other three.
!
    p(in1) = gfl(indx(3,in1))
    p(in2) = gfl(indx(3,in2))
    p(in3) = gfl(indx(3,in3))
    p(in4) = 0.5D+00 *(p(in1)+p(in2))
    p(in5) = 0.5D+00 *(p(in2)+p(in3))
    p(in6) = 0.5D+00 *(p(in3)+p(in1))

  end do

  return
end
subroutine l2norm(gfl,gflnrm,indx,nelem,neqnfl,node,np,xc,yc)

!*****************************************************************************80
!
!! L2NORM computes the "big" L2 norm of the velocity over the flow region.
!
!
!  Discussion:
!
!    A 13 point Gauss rule is used.
!
!    Note that this is the "BIG L2" norm, that is, the square root
!    of the integral of the square of the velocity over the flow region,
!    and NOT the "little l2" norm, which is simply the square root of the
!    sum of the squares of the coefficients.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!
!         GFL is the current solution vector, in which are stored
!         the finite element coefficients that define the velocity
!         and pressure functions, U, V and P.
!
!  GFLNRM Output, double recision GFLNRM.
!         GFLNRM is the approximate value of the square root of
!         the integral of the square of the velocity over the
!         flow domain.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of finite element equations used
!         to define the horizontal and vertical velocities and the
!         pressure.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM).
!
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  XC     Input, real ( kind = 8 ) XC(NP).
!
!         The X coordinates of the nodes.
!
!  YC     Input, real ( kind = 8 ) YC(NP).
!
!         The Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) nquad
!
  parameter (nquad = 13)
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) area
  real ( kind = 8 ) area2
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) eta
  real ( kind = 8 ) etaquad(nquad)
  real ( kind = 8 ) gfl(neqnfl)
  real ( kind = 8 ) gflnrm
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) in
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) jp
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) vmax
  real ( kind = 8 ) w
  real ( kind = 8 ) wquad(nquad)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsiquad(nquad)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq
!
  wquad(1)  = 0.175615257433204D+00
  wquad(2)  = 0.175615257433204D+00
  wquad(3)  = 0.175615257433204D+00

  wquad(4)  = 0.053347235608839D+00
  wquad(5)  = 0.053347235608839D+00
  wquad(6)  = 0.053347235608839D+00

  wquad(7)  = 0.077113760890257D+00
  wquad(8)  = 0.077113760890257D+00
  wquad(9)  = 0.077113760890257D+00
  wquad(10) = 0.077113760890257D+00
  wquad(11) = 0.077113760890257D+00
  wquad(12) = 0.077113760890257D+00

  wquad(13) = -0.149570044467670D+00

  do i = 1,nquad
    wquad(i) = wquad(i) / 2.0D+00
  end do

  xsiquad(1) =  0.260345966079038D+00
  etaquad(1) =  0.479308067841923D+00

  xsiquad(2) =  0.260345966079038D+00
  etaquad(2) =  0.260345966079038D+00

  xsiquad(3) =  0.479308067841923D+00
  etaquad(3) =  0.260345966079038D+00

  xsiquad(4) =  0.065130102902216D+00
  etaquad(4) =  0.869739794195568D+00

  xsiquad(5) =  0.065130102902216D+00
  etaquad(5) =  0.065130102902216D+00

  xsiquad(6) =  0.869739794195568D+00
  etaquad(6) =  0.065130102902216D+00

  xsiquad(7) =  0.048690315425316D+00
  etaquad(7) =  0.638444188569809D+00

  xsiquad(8) =  0.312865496004875D+00
  etaquad(8) =  0.638444188569809D+00

  xsiquad(9) =  0.048690315425316D+00
  etaquad(9) =  0.312865496004875D+00

  xsiquad(10) = 0.638444188569809D+00
  etaquad(10) = 0.312865496004875D+00

  xsiquad(11) = 0.312865496004875D+00
  etaquad(11) = 0.048690315425316D+00

  xsiquad(12) = 0.638444188569809D+00
  etaquad(12) = 0.048690315425316D+00

  xsiquad(13) = 1.0D+00 / 3.0D+00
  etaquad(13) = 1.0D+00 / 3.0D+00

  do i = 1,nquad
    xsiquad(i) = 1.0D+00 - xsiquad(i)
  end do

  gflnrm = 0.0D+00
  area2 = 0.0D+00
  vmax = 0.0D+00
!
!  Consider an element.
!
  do ielem = 1,nelem
!
!  Compute the area of the element.  For now, we assume that all
!  elements are triangles, and NOT isoparametric!
!
    ip1 = node(1,ielem)
    ip2 = node(2,ielem)
    ip3 = node(3,ielem)

    area = abs( &
         (yc(ip1)+yc(ip2))*(xc(ip2)-xc(ip1)) &
        +(yc(ip2)+yc(ip3))*(xc(ip3)-xc(ip2)) &
        +(yc(ip3)+yc(ip1))*(xc(ip1)-xc(ip3)) )
!
!  Evaluate the integrand at the quadrature points.
!
    do iquad = 1,nquad

      xsi = xsiquad(iquad)
      eta = etaquad(iquad)
      call xofxsi(eta,ielem,nelem,node,np,xq,xc,xsi,yq,yc)
!
!  Evaluate U, V and P at the IQUAD-th quadrature point by
!  finding the value of each of the 6 basis functions there,
!  and multiplying by their respective coefficients.
!
      u = 0.0D+00
      v = 0.0D+00
      do in = 1,6
        call qbf(ielem,in,w,dwdx,dwdy,nelem,node,np,xc,xq,yc,yq)
        ip = node(in,ielem)
        jp = indx(1,ip)
        u = u+gfl(jp)*w
        jp = indx(2,ip)
        v = v+gfl(jp)*w
      end do

      gflnrm = gflnrm+area*wquad(iquad)*(u**2+v**2)
      if ( vmax < u**2 + v**2 ) then
        vmax = u**2+v**2
      end if

      area2 = area2+area*wquad(iquad)
    end do
  end do

  gflnrm = sqrt(gflnrm)

  return
end
subroutine lbase ( ival, npol, pval, xpol, xval )

!*****************************************************************************80
!
!! LBASE evalualates the IVAL-th Lagrange basis polynomial.
!
!
!  Discussion:
!
!    The Lagrange interpolation basis polynomials are based
!    on the NPOL points XPOL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, the polynomial to evaluate.
!    IVAL should be between 1 and NPOL.
!
!    Input, integer ( kind = 4 ) NPOL, the number of points that define
!    the Lagrange polynomials.
!
!    Output, real ( kind = 8 ) PVAL, the value of the IVAL-th
!    Lagrange polynomial at the point XVAL.
!
!    Input, real ( kind = 8 ) XPOL(NPOL), the abscissas of the
!    Lagrange polynomials.  The entries in XPOL should be
!    distinct.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the
!    IVAL-th Lagrange polynomial is to be evaluated.
!
  implicit none
!
  integer ( kind = 4 ) npol
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) pval
  real ( kind = 8 ) xpol(npol)
  real ( kind = 8 ) xval
!
  pval = 1.0D+00
  do i = 1, npol
    if ( i /= ival ) then
      pval = pval * ( xval - xpol(i) ) / ( xpol(ival) - xpol(i) )
    end if
  end do

  return
end
function s_eqidb(strng1,strng2)

!*****************************************************************************80
!
!! S_EQIDB is a case insensitive comparison of two strings for
!  equality, ignoring blanks. 
!
!  Thus, LEQIDB('Nor Way','NORway') is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  STRNG1,
!  STRNG2 Input, character*(*) STRNG1, STRNG2, the strings to
!         compare.
!
!  LEQIDB Output, logical LEQIDB, the result of the comparison.
!
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  logical s_eqidb
  character s1
  character s2
  character ( len = * ) strng1
  character ( len = * ) strng2
!
  len1 = len(strng1)
  len2 = len(strng2)

  s_eqidb = .false.
 
  i1 = 0
  i2 = 0

10    continue
!
!  If we've matched all the nonblank characters in both strings,
!  then return with LEQIDB = .TRUE.
!
  if ( i1 == len1.and.i2.eq.len2) then
    s_eqidb = .true.
    return
  end if
!
!  Get S1, the next nonblank character in the first string.
!
20    continue

  i1 = i1+1
  if ( len1 < i1 ) then
    return
  end if
 
  if ( strng1(i1:i1) == ' ')go to 20
  s1 = strng1(i1:i1)
!
!  Get S2, the next nonblank character in the second string.
!
30    continue

  i2 = i2+1
  if ( len2 < i2 ) then
    return
  end if 
  if ( strng2(i2:i2) == ' ')go to 30
  s2 = strng2(i2:i2)
 
  if ( s1 /= s2)return

  go to 10

end     
function nbinom(m,n)

!*****************************************************************************80
!
!! NBINOM calculates a binomial coefficient.
!
!
!  Discussion:
!
!    The routine calculates the number of combinations of M things taken N
!    at a time.  NBINOM is ACM algorithm 160 translated to FORTRAN. 
!
!    The formula used is
!
!      NBINOM = M! / ( N! * (M-N)! )
!
!    This value is calculated in a way that tries to avoid overflow.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M.
!    M is the number of objects to choose from in the set,
!    or the row of Pascal's triangle in which the coefficient lies.
!    M should be zero or greater.
!
!  N      Input, integer ( kind = 4 ) N.
!         N is the number of objects selected from the set,
!         or the column of Pascal's triangle in which the coefficient
!         lies.  N should be 0 or greater, and no greater than M.
!
!  NBINOM Output, integer ( kind = 4 ) NBINOM.
!         NBINOM is the number of combinations of M things taken N
!         at a time.
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) n1
  integer ( kind = 4 ) nbinom
!
  if ( m < 0 ) then
    nbinom = 0
    return
  end if
 
  if ( n < 0 ) then
    nbinom = 0
    return
  end if
 
  if ( m < n ) then
    nbinom = 0
    return
  end if

  if ( n < m-n ) then
    n1 = m-n
  else
    n1 = n
  end if
 
  nbinom = 1

  do i = 1, m-n1
    nbinom = (nbinom*(n1+i)) / i
  end do
   
  return
end
subroutine nrmflo ( gfl, indx, neqnfl, np, resfl )

!*****************************************************************************80
!
!! NRMFLO returns norms of a flow solution or flow residual.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!         GFL is the current solution vector, in which are stored
!         the finite element coefficients that define the velocity
!         and pressure functions, U, V and P.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of finite element equations used
!         to define the horizontal and vertical velocities and the
!         pressure.
!
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
  implicit none
!
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) anrmf
  real ( kind = 8 ) anrmfp
  real ( kind = 8 ) anrmfu
  real ( kind = 8 ) anrmfv
  real ( kind = 8 ) anrmg
  real ( kind = 8 ) anrmp
  real ( kind = 8 ) anrmu
  real ( kind = 8 ) anrmv
  real ( kind = 8 ) enrmf
  real ( kind = 8 ) enrmfp
  real ( kind = 8 ) enrmfu
  real ( kind = 8 ) enrmfv
  real ( kind = 8 ) enrmg
  real ( kind = 8 ) enrmp
  real ( kind = 8 ) enrmu
  real ( kind = 8 ) enrmv
  real ( kind = 8 ) fp
  real ( kind = 8 ) fu
  real ( kind = 8 ) fv
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) inrmf
  integer ( kind = 4 ) inrmfp
  integer ( kind = 4 ) inrmfu
  integer ( kind = 4 ) inrmfv
  integer ( kind = 4 ) inrmg
  integer ( kind = 4 ) inrmp
  integer ( kind = 4 ) inrmu
  integer ( kind = 4 ) inrmv
  real ( kind = 8 ) p
  real ( kind = 8 ) resfl(neqnfl)
  real ( kind = 8 ) u
  real ( kind = 8 ) v
!
  anrmf = 0.0D+00
  anrmfp = 0.0D+00
  anrmfu = 0.0D+00
  anrmfv = 0.0D+00
  anrmg = 0.0D+00
  anrmp = 0.0D+00
  anrmu = 0.0D+00
  anrmv = 0.0D+00
  enrmfp = 0.0D+00
  enrmf = 0.0D+00
  enrmfu = 0.0D+00
  enrmfv = 0.0D+00
  enrmg = 0.0D+00
  enrmp = 0.0D+00
  enrmu = 0.0D+00
  enrmv = 0.0D+00
  inrmf = 1
  inrmfp = 1
  inrmfu = 1
  inrmfv = 1
  inrmg = 1
  inrmp = 1
  inrmu = 1
  inrmv = 1

  do i = 1,np

    u = gfl(indx(1,i))
    enrmu = enrmu+u**2
    enrmg = enrmg+u**2
    if ( anrmu < abs(u) ) then
      anrmu = abs(u)
      inrmu = i
    end if
    if ( anrmg < abs(u) ) then
      anrmg = abs(u)
      inrmg = i
    end if

    fu = resfl(indx(1,i))
    enrmf = enrmf+fu**2
    enrmfu = enrmfu+fu**2
    if ( anrmf < abs(fu) ) then
      anrmf = abs(fu)
      inrmf = i
    end if
    if ( anrmfu < abs(fu) ) then
      anrmfu = abs(fu)
      inrmfu = i
    end if

    v = gfl(indx(2,i))
    enrmv = enrmv+v**2
    enrmg = enrmg+v**2
    if ( anrmv < abs(v) ) then
      anrmv = abs(v)
      inrmv = i
    end if
    if ( anrmg <= abs(v) ) then
      anrmg = abs(v)
      inrmg = i
    end if

    fv = resfl(indx(2,i))
    enrmf = enrmf+fv**2
    enrmfv = enrmfv+fv**2
    if ( anrmf < abs(fv) ) then
      anrmf = abs(fv)
      inrmf = i
    end if
    if ( anrmfv < abs(fv) ) then
      anrmfv = abs(fv)
      inrmfv = i
    end if

    if ( 0 < indx(3,i) ) then

      p = gfl(indx(3,i))
      enrmp = enrmp+p**2
      enrmg = enrmg+p**2
      if ( anrmp <= abs(p) ) then
        inrmp = i
        anrmp = abs(p)
      end if
      if ( anrmg <= abs(p) ) then
        inrmg = i
        anrmg = abs(p)
      end if

      fp = resfl(indx(3,i))
      enrmf = enrmf+fp**2
      enrmfp = enrmfp+fp**2
      if ( anrmf < abs(fp) ) then
        inrmf = i
        anrmf = abs(fp)
      end if
      if ( anrmfp < abs(fp) ) then
        anrmfp = abs(fp)
        inrmfp = i
      end if

    end if

  end do

  enrmf = sqrt(enrmf)
  enrmfp = sqrt(enrmfp)
  enrmfu = sqrt(enrmfu)
  enrmfv = sqrt(enrmfv)
  enrmg = sqrt(enrmg)
  enrmp = sqrt(enrmp)
  enrmu = sqrt(enrmu)
  enrmv = sqrt(enrmv)
!
!  Print out results.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     MxNorm       Node    l2 Norm'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,i6,g14.6)' ) 'U  ', anrmu,  inrmu,  enrmu
  write ( *, '(a,g14.6,i6,g14.6)' ) 'V  ', anrmv,  inrmv,  enrmv
  write ( *, '(a,g14.6,i6,g14.6)' ) 'P  ', anrmp,  inrmp,  enrmp
  write ( *, '(a,g14.6,i6,g14.6)' ) 'GFL', anrmg,  inrmg,  enrmg
  write ( *, '(a,g14.6,i6,g14.6)' ) 'FU ', anrmfu, inrmfu, enrmfu
  write ( *, '(a,g14.6,i6,g14.6)' ) 'FV ', anrmfv, inrmfv, enrmfv
  write ( *, '(a,g14.6,i6,g14.6)' ) 'FP ', anrmfp, inrmfp, enrmfp
  write ( *, '(a,g14.6,i6,g14.6)' ) 'F  ', anrmf,  inrmf,  enrmf

  return
end
subroutine pcval ( nvec, xval, xvec, yval, yvec )

!*****************************************************************************80
!
!! PCVAL evaluates a piecewise constant function at a given point.
!
!
!  Discussion:
!
!    The piecewise constant function is specified as suggested by the
!    following graph:
!
!
!  Y(2)->              *---------*
!                      |         |
!  Y(1)->   *----------*         |
!                                |
!                                |
!  Y(3)->                        *---------......
!                             
!           ^          ^         ^
!           |          |         |
!         X(1)        X(2)      X(3)
!
!  Note that if XVAL falls to the left of XVEC(1), then YVAL = YVEC(1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise constant function.  NVEC must be at
!    least 1.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is to be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YVAL, the value of the piecewise
!    constant function at the point XVAL.
!
!    Input, real ( kind = 8 ) YVEC(NVEC), the value of the piecewise
!    constant function at each of the abscissas.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
  real ( kind = 8 ) yvec(nvec)
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1)) then
    yval = yvec(1)
    return
  else if ( xvec(nvec) <= xval ) then
    yval = yvec(nvec)
    return
  end if
!
!  Step 2: Find index I so that XVEC(I) <= XVAL < XVEC(I+1)
!
  do i = 1, nvec-1

    if ( xvec(i) <= xval .and. xval <= xvec(i+1) ) then
      yval = xvec(i)
      return
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PCVal - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ', xval
  stop
end
subroutine pldx ( nvec, xval, xvec, yder, yvec )

!*****************************************************************************80
!
!! PLDX evaluates the derivative of a piecewise linear function.
!
!
!  Discussion:
!
!    Note that if XVAL falls to the left of XVEC(1), then YDER = 0,
!    and similarly, if XDER is greater than XVEC(NVEC), YVAL = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise linear.  NVEC must be odd, and
!    at least 3.
!
!  XVAL   Input, real ( kind = 8 ) XVAL, the point at which the
!         derivative with respect to X is to be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!         function.  These should be distinct and in ascending order.
!
!  YDER   Output, real ( kind = 8 ) YDER, the value of the derivative of
!         the piecewise linear function with respect to X, at the point
!         XVAL.
!
!  YVEC   Input, real ( kind = 8 ) YVEC(NVEC), the value of the
!         piecewise linear function at each of the abscissas.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
  real ( kind = 8 ) yvec(nvec)
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1)) then
    yder = 0
    return
  else if ( xvec(nvec) <= xval ) then
    yder = 0
    return
  end if
!
!  Step 2: Find index I so that XVEC(I) <= XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i) <= xval.and.xval <= xvec(i+1)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLVal - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Evaluate the slope of the linear function at XVAL.
!
  i = ival

  yder = (yvec(i+1)-yvec(i)) / (xvec(i+1)-xvec(i))

  return
end
subroutine pldx1(ivec,nvec,xval,xvec,yder)

!*****************************************************************************80
!
!! PLDX1 evaluates the X derivative of the piecewise linear
!  polynomial which is 1 at the IVEC-th node and 0 at the others.
!
!  Note that if XVAL falls to the left of XVEC(1), then YDER = 0,
!  and similarly, if XVAL is greater than XVEC(NVEC), YDER = 0.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!    the partial derivative is desired.
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise linear.  NVEC must be odd, and
!    at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is to be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YDER, the value of the derivative of
!    the piecewise linear function at the point XVAL.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1)) then
    yder = 0.0D+00
    return
  else if ( xvec(nvec) <= xval ) then
    yder = 0.0D+00
    return
  end if
!
!  Step 2: Find index I so that XVEC(I) <= XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i) <= xval.and.xval <= xvec(i+1)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLDX1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Evaluate the slope of the IVEC-th linear function at XVAL.
!
  i = ival
  if ( ival == ivec) then
    yder = (0.0D+00 -1.0D+00 )/(xvec(ival+1)-xvec(ival))
  else if ( ival+1 == ivec) then
    yder = (1.0D+00 -0.0D+00 )/(xvec(ival+1)-xvec(ival))
  else
    yder = 0.0D+00
  end if

  return
end
subroutine pltopn ( disfil, igunit )

!*****************************************************************************80
!
!! PLTOPN opens the plotting file.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 30 ) DISFIL, the name of the file into 
!    which the graphics information will be stored.
!
!  IGUNIT Input/output, integer ( kind = 4 ) IGUNIT.
!
!         On input, if IGUNIT is zero, then the routine believes
!         that the graphics unit has not yet been opened.
!
!         If the FORTRAN unit has already been opened, then IGUNIT
!         should be nonzero, and the routine will know not to try
!         to open the file, since it is already open.
!
!         On output, IGUNIT is the FORTRAN unit used for writing data
!         to the plotfile FILEG.
!
  implicit none
!
  character ( len = 30 ) disfil
  integer ( kind = 4 ) igunit
!
!  If IGUNIT is not zero, then the graphics unit has already
!  been opened.
!
  if ( igunit == 0) then

    write ( *, * ) ' '
    write ( *, * ) 'PltOpn - Note:'
    write ( *, * ) '  Opening the DISPLAY plot file '//disfil
    write ( *, * ) ' '
!
!  Delete any old copy of the file.
!
    igunit = 11
    open(unit = igunit,file=disfil,status='unknown', &
         form = 'formatted',access='sequential',err=10)

    return
!
!  Write a warning if the plot file could not be opened.
!
10      continue

    write ( *, * ) ' '
    write ( *, * ) 'PltOpn - Warning!'
    write ( *, * ) '  The plot file could not be opened.'
    igunit = 0
!
!  Else plotfile is already open.
!
  else
    write ( *, * ) ' '
    write ( *, * ) 'PltOpn - Note'
    write ( *, * ) '  The plot file is already open.'
    write ( *, * ) '  New information will be appended to it.'
  end if

  return
end
subroutine plval ( nvec, xval, xvec, yval, yvec )

!*****************************************************************************80
!
!! PLVAL evaluates a piecewise linear function at a given point.
!
!
!  Note that if XVAL falls to the left of XVEC(1), then YVAL = YVEC(1),
!  and similarly, if XVAL is greater than XVEC(NVEC), YVAL = YVEC(NVEC).
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  NVEC   Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!         that define the piecewise linear.  NVEC must be at least 1.
!
!  XVAL   Input, real ( kind = 8 ) XVAL, the point at which the function
!         is to be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!         function.  These should be distinct and in ascending order.
!
!  YVAL   Output, real ( kind = 8 ) YVAL, the value of the piecewise
!         linear function at the point XVAL.
!
!  YVEC   Input, real ( kind = 8 ) YVEC(NVEC), the value of the piecewise
!         function at each of the abscissas.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
  real ( kind = 8 ) yvec(nvec)
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1)) then
    yval = yvec(1)
    return
  else if ( xvec(nvec) <= xval ) then
    yval = yvec(nvec)
    return
  end if
!
!  Step 2: Find index I so that XVEC(I) <= XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i) <= xval.and.xval <= xvec(i+1)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLVal - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Evaluate the linear function at XVAL.
!
  i = ival

  if ( xval == xvec(i+1)) then
    yval = yvec(i+1)
  else if ( xval == xvec(i)) then
    yval = yvec(i)
  else
    yval = ( yvec(i)*(xvec(i+1)-xval) &
           +yvec(i+1)*(xval-xvec(i)) ) / (xvec(i+1)-xvec(i))
  end if

  return
end
subroutine plval1 ( ivec, nvec, xval, xvec, yval )

!*****************************************************************************80
!
!! PLVAL1 evaluates a piecewise linear basis polynomial.
!
!
!  Discussion:
!
!    The piecewise linear basis polynomial is 1
!    at node IVEC and 0 at the other nodes.
!
!    Note that if XVAL falls to the left of XVEC(1), then YVAL = 0,
!    and similarly, if XVAL is greater than XVEC(NVEC), YVAL = 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!    the partial derivative is desired.
!
!    Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!    that define the piecewise linear.  NVEC must be odd, and
!    at least 3.
!
!    Input, real ( kind = 8 ) XVAL, the point at which the function
!    is to be evaluated.
!
!    Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!    function.  These should be distinct and in ascending order.
!
!    Output, real ( kind = 8 ) YDER, the value of the derivative of
!    the piecewise linear function at the point XVAL.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
!
!  Step 1: Check if XVAL lies outside the intervals.
!
  if ( xval <= xvec(1)) then
    yval = 0.0D+00
    return
  else if ( xvec(nvec) <= xval ) then
    yval = 0.0D+00
    return
  end if
!
!  Step 2: Find index I so that XVEC(I) <= XVAL < XVEC(I+1)
!
  do i = 1,nvec-1

    if ( xvec(i) <= xval.and.xval <= xvec(i+1)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PLVAL1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 3: Determine the index of the left endpoint of the least and
!  greatest intervals that IVEC can affect.
!
  i = ival
  if ( ival == ivec) then
    if ( xval == xvec(ival)) then
      yval = 1.0D+00
    else
      yval = (xvec(ival+1)-xval) / (xvec(ival+1)-xvec(ival))
    end if
  else if ( ival+1 == ivec) then
    if ( xval == xvec(ival+1) ) then
      yval = 1.0D+00
    else
      yval = (xval-xvec(ival)) / (xvec(ival+1)-xvec(ival))
    end if
  else
    yval = 0.0D+00
  end if

  return
end
subroutine pqdx ( nvec, xval, xvec, yder, yvec )

!*****************************************************************************80
!
!! PQDX evaluates the derivative of a piecewise quadratic function with
!  respect to its argument at a given point.
!
!  Note that if XDER falls to the left of XVEC(1), then YVAL = 0,
!  and similarly, if XVAL is greater than XVEC(NVEC), YDER = 0.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  NVEC   Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!         that define the piecewise quadratic.  NVEC must be odd, and
!         at least 3.
!
!  XVAL   Input, real ( kind = 8 ) XVAL, the point at which the
!         derivative with respect to X is to be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!         function.  These should be distinct and in ascending order.
!
!  YDER   Output, real ( kind = 8 ) YDER, the value of the derivative
!         of the piecewise  quadratic function with respect to X,
!         at the point XVAL.
!
!  YVEC   Input, real ( kind = 8 ) YVEC(NVEC), the value of the piecewise
!         quadratic function at each of the abscissas.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
  real ( kind = 8 ) yvec(nvec)
!
!  Step 0: Check data.
!
  if ( nvec < 3) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX - Fatal error.'
    write ( *, * ) '  NVEC is ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2) /= 1) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I) <= XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1)) then
    yder = yvec(1)
    return
  else if ( xvec(nvec) <= xval ) then
    yder = yvec(nvec)
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval.and.xval <= xvec(i+2)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQDX - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 2: Evaluate the derivative of the quadratic function at XVAL.
!
  i = ival

  yder = yvec(i)*(2*xval-xvec(i+1)-xvec(i+2)) &
    /((xvec(i)-xvec(i+1))*(xvec(i)-xvec(i+2))) &
         +yvec(i+1)*(2*xval-xvec(i)-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2))) &
         +yvec(i+1)*(2*xval-xvec(i)-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2)))

  return
end
subroutine pqdx1 ( ivec, nvec, xval, xvec, yder )

!*****************************************************************************80
!
!! PQDX1 evaluates the X derivative of the piecewise quadratic
!  polynomial which is 1 at the IVEC-th node and 0 at the others.
!
!  Note that if XVAL falls to the left of XVEC(1), then YDER = 0,
!  and similarly, if XVAL is greater than XVEC(NVEC), YDER = 0.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IVEC   Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!         the partial derivative is desired.
!
!  NVEC   Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!         that define the piecewise quadratic.  NVEC must be odd, and
!         at least 3.
!
!  XVAL   Input, real ( kind = 8 ) XVAL, the point at which the function
!         be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!         function.  These should be distinct and in ascending order.
!
!  YDER   Output, real ( kind = 8 ) YDER, the value of the derivative of
!         the piecewise quadratic function at the point XVAL.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yder
!
!  Step 0: Check data.
!
  if ( nvec < 3) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX1 - Fatal error!'
    write ( *, * ) '  NVEC = ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2) /= 1) then
    write ( *, * ) ' '
    write ( *, * ) 'PQDX1 - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I) <= XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1)) then
    yder = 0
    return
  else if ( xvec(nvec) <= xval ) then
    yder = 0
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval.and.xval <= xvec(i+2)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQDX1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 2: Determine the index of the left endpoint of the least and
!  greatest intervals that IVEC can affect.
!
  if ( mod(ivec,2) == 0) then
    ilo = ivec-1
    ihi = ivec-1
  else
    ilo = max(ivec-2,1)
    ihi = ivec
  end if
!
!  Step 3: If XVAL is outside of the intervals that IVEC can affect,
!  the derivative is zero.
!
  if ( ival < ilo .or. ihi < ival ) then
    yder = 0
    return
  end if
!
!  Step 3: Evaluate the derivative of the quadratic function at XVAL.
!
  i = ival

  if ( ivec == ival ) then
    yder = ( 2.0D+00 * xval - xvec(i+1) - xvec(i+2) ) &
           / ( ( xvec(i) - xvec(i+1) ) * ( xvec(i) - xvec(i+2) ) )
  else if ( ivec == ival+1) then
     yder = ( 2.0D+00 * xval-xvec(i)-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2)))
  else if ( ivec == ival+2) then
      yder = ( 2.0D+00 * xval-xvec(i)-xvec(i+1)) &
        /((xvec(i+2)-xvec(i))*(xvec(i+2)-xvec(i+1)))
  else
    write ( *, * ) ' '
    write ( *, * ) 'PQDX1 - Fatal error!'
    write ( *, * ) '  IVEC = ',ivec
    write ( *, * ) '  IVAL = ',ival
  end if

  return
end
subroutine pqval ( nvec, xval, xvec, yval, yvec )

!*****************************************************************************80
!
!! PQVAL evaluates a piecewise quadratic function at a given point.
!
!
!  Discussion:
!
!    The piecewise quadratic is defined by NVEC values, where NVEC
!    is odd, and at least 3.  The function is defined by specifying
!    a list of nodes XVEC(I), and specifying its value YVEC(I) at each
!    node. 
!
!    The function will be a quadratic polynomial over each of
!    (NVEC-1)/2 intervals that are made up a set of three consecutive
!    nodes, with the first one odd.  Thus, XVEC(1), XVEC(2) and XVEC(3)
!    lie in the first interval.
!
!    At the odd nodes, the quadratic that defines the function may
!    change, but the function remains continuous there, though not
!    differentiable.
!
!    Note that if XVAL falls to the left of XVEC(1), then YVAL = YVEC(1),
!    and similarly, if XVAL is greater than XVEC(NVEC), YVAL = YVEC(NVEC).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  NVEC   Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!         that define the piecewise quadratic. 
!
!         NVEC must be odd, and at least 3.
!
!  XVAL   Input, real ( kind = 8 ) XVAL, the point at which the function
!         is be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!         function.  These should be distinct and in ascending order.
!
!  YVAL   Output, real ( kind = 8 ) YVAL, the value of the piecewise
!         quadratic function at the point XVAL.
!
!  YVEC   Input, real ( kind = 8 ) YVEC(NVEC), the value of the
!         piecewise quadratic function at each of the abscissas.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ival
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
  real ( kind = 8 ) yvec(nvec)
!
!  Step 0: Check data.
!
  if ( nvec < 3) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVal - Fatal error!'
    write ( *, * ) '  Value of NVEC = ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2) /= 1) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVal - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I) <= XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1)) then
    yval = yvec(1)
    return
  else if ( xvec(nvec) <= xval ) then
    yval = yvec(nvec)
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval.and.xval <= xvec(i+2)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQVal - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  write ( *, * ) '  There are ',nvec,' nodes.'
  write ( *, * ) '  First node is at ',xvec(1)
  write ( *, * ) '  Last node is at  ',xvec(nvec)

  do i = 1,nvec
    write ( *, * ) xvec(i)
  end do
  stop

10    continue
!
!  Step 2: Evaluate the quadratic function at XVAL.
!
  i = ival

  yval = yvec(i)*(xval-xvec(i+1)) * (xval-xvec(i+2)) &
         /((xvec(i)-xvec(i+1))*(xvec(i)-xvec(i+2))) &
         +yvec(i+1)*(xval-xvec(i)) * (xval-xvec(i+2)) &
         /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2))) &
         +yvec(i+2)*(xval-xvec(i)) * (xval-xvec(i+1)) &
         /((xvec(i+2)-xvec(i))*(xvec(i+2)-xvec(i+1)))

  return
end
subroutine pqval1 ( ivec, nvec, xval, xvec, yval )

!*****************************************************************************80
!
!! PQVAL1 evaluates the piecewise quadratic polynomial which is 1
!  at node IVEC and 0 at the other nodes.
!
!  Note that if XVAL falls to the left of XVEC(1), then YVAL = 0,
!  and similarly, if XVAL is greater than XVEC(NVEC), YVAL = 0.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IVEC   Input, integer ( kind = 4 ) IVEC, the coefficient with respect to which
!         the partial derivative is desired.
!
!  NVEC   Input, integer ( kind = 4 ) NVEC, the number of abscissas and coefficients
!         that define the piecewise quadratic.  NVEC must be odd, and
!         at least 3.
!
!  XVAL   Input, real ( kind = 8 ) XVAL, the point at which the function
!         is to be evaluated.
!
!  XVEC   Input, real ( kind = 8 ) XVEC(NVEC), the abscissas of the
!         function.  These should be distinct and in ascending order.
!
!  YDER   Output, real ( kind = 8 ) YDER, the value of the derivative of
!         the piecewise quadratic function at the point XVAL.
!
  implicit none
!
  integer ( kind = 4 ) nvec
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) ivec
  real ( kind = 8 ) xval
  real ( kind = 8 ) xvec(nvec)
  real ( kind = 8 ) yval
!
!  Step 0: Check data.
!
  if ( nvec < 3) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVal1 - Fatal error!'
    write ( *, * ) '  Value of NVEC is ',nvec
    write ( *, * ) '  but NVEC must be at least 3.'
    stop
  end if

  if ( mod(nvec,2) /= 1) then
    write ( *, * ) ' '
    write ( *, * ) 'PQVal1 - Fatal error!'
    write ( *, * ) '  Even value of NVEC = ',nvec
    stop
  end if
!
!  Step 1: Find odd index I so that XVEC(I) <= XVAL < XVEC(I+2)
!
  if ( xval <= xvec(1)) then
    yval = 0.0D+00
    return
  else if ( xvec(nvec) <= xval ) then
    yval = 0.0D+00
    return
  end if

  do i = 1,nvec-2,2

    if ( xvec(i) <= xval.and.xval <= xvec(i+2)) then
      ival = i
      go to 10
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'PQVal1 - Fatal error!'
  write ( *, * ) '  Could not bracket XVAL = ',xval
  stop

10    continue
!
!  Step 2: Determine the index of the left endpoint of the least and
!  greatest intervals that IVEC can affect.
!
  if ( mod(ivec,2) == 0) then
    ilo = ivec-1
    ihi = ivec-1
  else
    ilo = max(ivec-2,1)
    ihi = ivec
  end if
!
!  Step 3: If XVAL is outside of the intervals that IVEC can affect,
!  the value is zero.
!
  if ( ival < ilo .or. ihi < ival ) then
    yval = 0
    return
  end if
!
!  Step 3: Evaluate the quadratic function at XVAL.
!
  i = ival

  if ( ivec == ival) then
    yval = (xval-xvec(i+1)) * (xval-xvec(i+2)) &
         /((xvec(i)-xvec(i+1))*(xvec(i)-xvec(i+2)))
  else if ( ivec == ival+1) then
     yval = (xval-xvec(i)) * (xval-xvec(i+2)) &
        /((xvec(i+1)-xvec(i))*(xvec(i+1)-xvec(i+2)))
  else if ( ivec == ival+2) then
      yval = (xval-xvec(i)) * (xval-xvec(i+1)) &
        /((xvec(i+2)-xvec(i))*(xvec(i+2)-xvec(i+1)))
  else
    write ( *, * ) ' '
    write ( *, * ) 'PQVal1 - Fatal error!'
    write ( *, * ) '  IVEC = ',ivec
    write ( *, * ) '  IVAL = ',ival
  end if

  return
end
subroutine prbmat ( afl, ihi, ilo, jhi, jlo, ldafl, neqnfl, nlband )

!*****************************************************************************80
!
!! PRBMAT prints the nonzeroes in a submatrix of a band matrix.
!
!
!  Discussion:
!
!    The submatrix is the rectangular region including rows ILO to IHI, 
!    columns JLO to JHI.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  AFL    Input, real ( kind = 8 ) AFL(LDAFL,MAXNFL).
!         If Newton iteration is being carried out, then AFL contains the
!         Jacobian matrix for the full system.  If Picard iteration is
!         being carried out, AFL contains the Picard matrix.
!         AFL is stored in LINPACK general band storage mode, with
!         dimension 3*NBANDL+1 by NEQNFL.
!
!  IHI,
!  ILO,
!  JHI,
!  JLO    Input, integer ( kind = 4 ) IHI, ILO, JHI, JLO.
!         PRMAT is to print all nonzero entries in rows ILO through IHI,
!         and columns JLO through JHI, of the matrix AFL.
!
!  LDAFL  Input, integer ( kind = 4 ) LDAFL.
!         LDAFL is the first dimension of the matrix AFL as declared in
!         the main program.  LDAFL must be at least 3*NLBAND+1.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL.
!         NEQNFL is the number of equations (and coefficients) in the full
!         finite element system.
!
!  NLBAND Input, integer ( kind = 4 ) NLBAND, the lower bandwidth of the matrix AFL. 
!         The zero structure of AFL is assumed to be symmetric, and so
!         NLBAND is also the upper bandwidth of AFL. 
!
  implicit none
!
  integer ( kind = 4 ), parameter :: incx = 5
!
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) neqnfl
!
  real ( kind = 8 ) afl(ldafl,neqnfl)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) nlband
!
  do j2lo = jlo,jhi,incx

    j2hi = j2lo+incx-1
    j2hi = min(j2hi,neqnfl)
    j2hi = min(j2hi,jhi)

    inc = j2hi+1-j2lo

    write ( *, * ) ' '
    do j = j2lo,j2hi
      j2 = j+1-j2lo
      write(ctemp(j2),'(i7,7x)')j
    end do
    write(*,'(''Columns '',5a14)')(ctemp(j2),j2 = 1,inc)
!       write ( *, * ) 'Columns ',j2lo,' to ',j2hi
    write ( *, * ) '  Row'
    write ( *, * ) ' '

    i2lo = ilo
    i2lo = max(ilo,1)
    i2lo = max(j2lo-nlband,i2lo)

    i2hi = ihi
    i2hi = min(ihi,neqnfl)
    i2hi = min(j2hi+nlband,i2hi)

    do i = i2lo,i2hi

      do j2 = 1,inc

        j = j2lo-1+j2

        if ( i-j <= nlband.and.j-i <= nlband) then
          write(ctemp(j2),'(g14.6)')afl(i-j+2*nlband+1,j)
          if ( afl(i-j+2*nlband+1,j) == 0.0D+00 )ctemp(j2)='    0.0'
        else
          ctemp(j2) = '              '
        end if

      end do

      write(*,'(i5,1x,5a14)')i,(ctemp(j2),j2 = 1,inc)

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine prdat(disfil,drey,epsdif,gridx,gridy,hx,hy,ibs,ibump,ifs,ijac, &
  iopt,maxnew,maxopt,maxsim,nbcrb,ncofrb,nelem,nferb,neqnfl,np,npar,nparb, &
  nparf,ntay,nx,ny,region,reytay,tecfil,tolnew,tolopt,tolsim,wateb,watep, &
  wateu,watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)

!*****************************************************************************80
!
!! PRDAT prints the problem information.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) npar
!
  character ( len = 30 ) disfil
  real ( kind = 8 ) drey
  real ( kind = 8 ) epsdif
  character ( len = 20 ) gridx
  character ( len = 20 ) gridy
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) ijac
  integer ( kind = 4 ) iopt(npar)
  integer ( kind = 4 ) maxnew
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) maxsim
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nferb
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) ntay
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  character ( len = 20 ) region
  real ( kind = 8 ) reytay
  character ( len = 30 ) tecfil
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  character ( len = 6 ) type
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  character ( len = 3 ) yesno
  real ( kind = 8 ) yrange
!
  write ( *, * ) ' '
  write ( *, * ) '  DISPLAY graphics file is DISFIL =              ', trim ( disfil )
  write ( *, * ) '  REYNLD increment for finite differences DREY = ',drey
  write ( *, * ) '  Finite difference perturbation EPSDIF =        ',epsdif
  write ( *, * ) '  X grid generation option GRIDX = '//gridx
  write ( *, * ) '  Y grid generation option GRIDY = '//gridy
  write ( *, * ) '  X spacing, HX =   ',hx
  write ( *, * ) '  Y spacing, HY =   ',hy
  write ( *, * ) '  Bump piecewise polynomial order IBS = ',ibs
  write ( *, * ) '  Bump option IBUMP =                   ',ibump
  write ( *, * ) '  Flow piecewise polynomial order IFS = ',ifs
  write ( *, * ) '  Jacobian option IJAC =                ',ijac
  write ( *, * ) ' '
  write ( *, * ) '  Variable  Type  Free to Vary?'
  write ( *, * ) ' '
  do i = 1,npar
    if ( i <= nparf) then
      type = 'Inflow'
    else if ( i <= nparf+nparb) then
      type = 'Shape'
    else
      type = 'Reynld'
    end if
    if ( iopt(i) == 0) then
      yesno = 'No'
    else
      yesno = 'Yes'
    end if
    write(*,'(6x,i5,2x,a6,2x,a3)')i,type,yesno
  end do
  write ( *, * ) ' '
  write ( *, * ) '  Maximum Newton iterations MAXNEW =    ',maxnew
  write ( *, * ) '  Maximum optimization steps MAXOPT =   ',maxopt
  write ( *, * ) '  Maximum Newton iterations MAXSIM =    ',maxsim
  write ( *, * ) '  # of RB boundary conditions NBCRB =   ',nbcrb
  write ( *, * ) '  Number of reduced equations, NCOFRB = ',ncofrb
  write ( *, * ) '  Number of elements, NELEM =           ',nelem
  write ( *, * ) '  Number of full equations, NEQNFL =    ',neqnfl
  write ( *, * ) '  # of FE reduced basis cofs, NFERB =   ',nferb
  write ( *, * ) '  Number of nodes, NP =                 ',np
  write ( *, * ) '  Number of parameters NPAR =           ',npar
  write ( *, * ) '  Number of inflow parameters NPARF =   ',nparf
  write ( *, * ) '  Number of Taylor vectors NTAY =       ',ntay
  write ( *, * ) '  Number of bump parameters NPARB =     ',nparb
  write ( *, * ) '  Number of X elements, NX =            ',nx
  write ( *, * ) '  Number of Y elements, NY =            ',ny
  write ( *, * ) '  The flow region is REGION = ', trim ( region )
  write ( *, * ) '  REYNLD value for Taylor, REYTAY =     ',reytay
  write ( *, * ) '  TECPLOT graphics file is TECFIL =              ', trim ( tecfil )
  write ( *, * ) '  Newton convergence tolerance TOLNEW = ',tolnew
  write ( *, * ) '  Optimization tolerance TOLOPT =       ',tolopt
  write ( *, * ) '  Picard convergence tolerance TOLSIM = ',tolsim
  write ( *, * ) '  Bump control cost,   WATEB =          ',wateb
  write ( *, * ) '  Pressure discrepancy, WATEP =         ',watep
  write ( *, * ) '  U discrepancy, WATEU =                ',wateu
  write ( *, * ) '  V discrepancy, WATEV =                ',watev
  write ( *, * ) '  Left X of bump, XBL =                 ',xbl
  write ( *, * ) '  Right X of bump, XBR =                ',xbr
  write ( *, * ) '  Flow profile measured at XPROF =      ',xprof
  write ( *, * ) '  X range, XRANGE =                     ',xrange
  write ( *, * ) '  Left Y of bump, YBL =                 ',ybl
  write ( *, * ) '  Right Y of bump, YBR =                ',ybr
  write ( *, * ) '  Y range, YRANGE =                     ',yrange

  return
end
subroutine prdmat(a,ihi,ilo,jhi,jlo,mhi,mlo,nhi,nlo)

!*****************************************************************************80
!
!! PRDMAT prints out a portion of a dense matrix.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real A(MLO:MHI,NLO:NHI), the matrix to be printed.
!
!  IHI,
!  ILO    Input, integer ( kind = 4 ) IHI, ILO.
!         ILO is the first and IHI the last row to print.
!
!  JHI,
!  JLO    Input, integer ( kind = 4 ) JHI, JLO.
!         JLO is the first, and JHI the last column to print.
!
!  MHI,
!  MLO    Input, integer ( kind = 4 ) MHI, MLO.
!         The rows of A go from MLO to MHI.
!
!  NHI,
!  NLO    Input, integer ( kind = 4 ) NHI, NLO.
!         The columns of A go from NLO to NHI.
!
  implicit none
!
  integer ( kind = 4 ) incx
  parameter (incx = 5)
!
  integer ( kind = 4 ) mhi
  integer ( kind = 4 ) mlo
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) nlo
!
  real ( kind = 8 ) a(mlo:mhi,nlo:nhi)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
!
  write ( *, * ) ' '

  do j2lo = jlo,jhi,incx

    j2hi = j2lo+incx-1
    if ( nhi < j2hi ) then
      j2hi = nhi
    end if

    if ( jhi < j2hi ) then
      j2hi = jhi
    end if

    inc = j2hi+1-j2lo

    write ( *, * ) ' '
    do j = j2lo,j2hi
      j2 = j+1-j2lo
      write(ctemp(j2),'(i7,7x)')j
    end do
    write(*,'(''Columns '',5a14)')(ctemp(j2),j2 = 1,inc)
    write ( *, * ) '  Row'
    write ( *, * ) ' '

    i2lo = max(ilo,mlo)
    i2hi = min(ihi,mhi)

    do i = i2lo,i2hi

      do j2 = 1,inc

        j = j2lo-1+j2
       
        write(ctemp(j2),'(g14.6)')a(i,j)
        if ( a(i,j) == 0.0D+00 )ctemp(j2)='    0.0'

      end do

      write(*,'(i5,1x,5a14)')i,(ctemp(j),j = 1,inc)

    end do

  end do

  write ( *, * ) ' '

  return
end
subroutine prelem ( ihi, ilo, nelem, node, np, xc, yc )

!*****************************************************************************80
!
!! PRELEM prints out data about one or more elements.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IHI, ILO.
!    ILO is the first element of interest, and IHI the last.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM.
!         NELEM is the numberof elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM), contains the numbers
!         of the nodes that make up each element.  Element number
!         I is associated with nodes NODE(1,I) through NODE(6,I).
!
!  NP     Input, integer ( kind = 4 ) NP.
!
!         NP is the number of nodes.  NP = (2*NX-1)*(2*NY-1).
!
!  XC,
!  YC     Input, real ( kind = 8 ) XC(NP), YC(NP).
!         XC and YC are the X and Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi2
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)
!
  if ( ilo < 1) then
    ilo2 = 1
  else
    ilo2 = ilo
  end if

  if ( nelem < ihi ) then
    ihi2 = nelem
  else
    ihi2 = ihi
  end if

  do ielem = ilo2,ihi2
    write ( *, * ) ' '
    write ( *, * ) 'Element IELEM = ',ielem
    write ( *, * ) ' '
    write ( *, * ) '  Nodes:'
    write ( *, * ) ' '
    do i = 1,6
      ip = node(i,ielem)
      write ( *, * ) i,ip,xc(ip),yc(ip)
    end do
  end do

  return
end
subroutine prfxfln ( neqnfl, resfl )

!*****************************************************************************80
!
!! PRFXFLN prints out the norm of a full residual.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NEQNFL, the number of equations (and coefficients) 
!    in the full finite element system.
!
!    Input, real ( kind = 8 ) RESFL(NEQNFL).
!    RESFL contains the residual in the full basis equations.
!
  implicit none
!
  integer ( kind = 4 ) neqnfl
!
  real ( kind = 8 ) anrmr
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) enrmr
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) itemp
  real ( kind = 8 ) resfl(neqnfl)
!
  itemp = idamax(neqnfl,resfl,1)
  anrmr = abs(resfl(itemp))
  enrmr = dnrm2(neqnfl,resfl,1)

  write ( *, * ) ' '
  write ( *, * ) '          MxNorm      l2 Norm'
  write ( *, * ) ' '
  write(*,'(''Fx(GFL) '',2g14.6)')anrmr,enrmr

  return
end
subroutine prgrb ( grb, ncofrb )

!*****************************************************************************80
!
!! PRGRB prints out the reduced basis solution.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) GRB(NCOFRB), coefficients for the
!    reduced system.
!
!    Input, integer ( kind = 4 ) NCOFRB, the number of coefficients for the
!    reduced system.
!
  implicit none
!
  integer ( kind = 4 ) ncofrb
!
  real ( kind = 8 ) anrmg
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) enrmg
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) itemp
!
  if ( ncofrb <= 0) then
    write ( *, * ) ' '
    write ( *, * ) 'PrGRB - Fatal error.'
    write ( *, * ) '  Input value of NCOFRB = ',ncofrb
    stop
  end if

  write ( *, * ) ' '
  write ( *, * ) 'PrGRB - The reduced basis coefficients:'
  write ( *, * ) ' '
  do i = 1,ncofrb
    write(*,'(i6,g14.6)')i,grb(i)
  end do

  itemp = idamax(ncofrb,grb,1)
  anrmg = abs(grb(itemp))
  enrmg = dnrm2(ncofrb,grb,1)

  write ( *, * ) ' '
  write ( *, * ) '          MxNorm      l2 Norm'
  write ( *, * ) ' '
  write(*,'(''GRB      '',2g14.6)')anrmg,enrmg

  return
end
subroutine prindx ( ihi, ilo, indx, np, xc, yc )

!*****************************************************************************80
!
!! PRINDX prints out the integer variables that define the
!  relationships between the nodes and elements.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IHI,
!  ILO    Input, integer ( kind = 4 ) IHI, ILO.
!         ILO is the first, and IHI the last node at which the
!         information is desired.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP).
!
!         INDX contains, for each node I, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry G(K).
!
!         If INDX(I,J) is positive, then that means that a degree of
!         freedom for variable J (U, V or P) is associated with node
!         I, and an equation will be generated to determine its value.
!
!         If INDX(I,J) is zero, then that means the the value of variabl
!         J (U, V or P) has been specified at node I.  No equation is
!         generated to determine its value.
!
!  NP     Input, integer ( kind = 4 ) NP.
!
!         NP is the number of nodes.  NP = (2*NX-1)*(2*NY-1).
!
!  XC,
!  YC     Input, real ( kind = 8 ) XC(NP), YC(NP).
!         XC and YC are the X and Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) np
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) indx(3,np)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)
!
  write ( *, * ) ' '
  write ( *, * ) 'PrIndx:'
  write ( *, * ) ' '
  write ( *, * ) ' Node    X             Y              U     V     P'
  write ( *, * ) ' '
  do i = max(ilo,1),min(ihi,np)
    if ( indx(3,i) /= 0) then
      write(*,'(i6,2g14.6,3i6)')i,xc(i),yc(i),indx(1,i),indx(2,i),indx(3,i)
    else
      write(*,'(i6,2g14.6,2i6)')i,xc(i),yc(i),indx(1,i),indx(2,i)
    end if
  end do

  return
end
subroutine prmatfl(a,eqn,ihi,ilo,indx,jhi,jlo,maxnfl,ncol,neqnfl,np)

!*****************************************************************************80
!
!! PRMATFL a matrix A associated with a full flow problem.
!
!
!  Discussion:
!
!    PRMATFL prints out a range of rows and columns of a dense matrix,
!    whose rows are indirectly indexed by node number, and whose
!    columns are indexed in the usual way.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  A      Input, real ( kind = 8 ) A(MAXNFL,NCOFRB).
!         A is the matrix whose entries are to be printed.
!
!  EQN    Input, character ( len = 2 ) EQN(MAXNFL).
!         EQN records the "type" of each equation that will be generated, and
!         which is associated with an unknown.  Note that most boundary
!         conditions do not result in an equation.  The current values are:
!
!         'U'  The horizontal momentum equation.
!         'UB' The condition U = 0 applied at a node on the bump.
!         'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!         'UW' The condition U = 0 applied at a node on a fixed wall.
!         'U0' A dummy value of U = 0 should be set.
!
!         'V'  The vertical momentum equation.
!         'VB' The condition V = 0 applied at a node on the bump.
!         'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!         'VW' The condition V = 0 applied at a node on a fixed wall.
!         'V0' A dummy value of V = 0 should be set.
!
!         'P'  The continuity equation.
!         'PB' The condition P = 0 applied at (XMAX,YMAX).
!         'P0' A dummy value of P = 0 should be set.
!
!  IHI,
!  ILO    Input, integer ( kind = 4 ) IHI, ILO.
!         ILO is the first node, and IHI the last node, for which the
!         data should be printed.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K).
! 
!         If INDX(I,J) is positive, then that means that a degree of
!         freedom for variable I (U, V or P) is associated with node
!         J, and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J.
!
!  JHI,
!  JLO    Input, integer ( kind = 4 ) JHI, JLO.
!         JLO is the first, and JHI the last column of A to print.
!
!  MAXNFL Input, integer ( kind = 4 ) MAXNFL.
!         MAXNFL is the maximum number of equations in the full system.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL.
!         NEQNFL is the number of equations (and coefficients) in the full
!         finite element system.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB.
!         NCOFRB is the number of sensitivities.
!
!  NP     Input, integer ( kind = 4 ) NP.
!         NP is the number of nodes.
!
  implicit none
!
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncol
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) a(maxnfl,ncol)
  character ( len = 2 ) eqn(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi2
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo2
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jhi2
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) jlo2
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lhi
  integer ( kind = 4 ) llo
  integer ( kind = 4 ) ncols
!
  if ( neqnfl <= 0) then
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Fatal error!'
    write ( *, * ) '  NEQNFL = ',neqnfl
    stop
  end if

  if ( ncol <= 0) then
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Fatal error!'
    write ( *, * ) '  NCOL = ',ncol
    stop
  end if

  if ( ilo < 1) then
    ilo2 = 1
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Warning!'
    write ( *, * ) '  Input ILO = ',ilo
    write ( *, * ) '  Reset to ILO2 = ',ilo2
  else
    ilo2 = ilo
  end if

  if ( np < ihi ) then
    ihi2 = np
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Warning!'
    write ( *, * ) '  Input IHI = ',ihi
    write ( *, * ) '  Reset to IHI2 = ',ihi2
  else
    ihi2 = ihi
  end if

  if ( ihi2 < ilo2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Warning:'
    write ( *, * ) '  Input ILO =      ',ilo, ' IHI= ',ihi
    write ( *, * ) '  Effective ILO2 = ',ilo2,' IHI2=',ihi2
    return
  end if

  if ( jlo < 1) then
    jlo2 = 1
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Warning!'
    write ( *, * ) '  Input value of JLO was ',jlo
    write ( *, * ) '  Reset to JLO2 = ',jlo2
  else
    jlo2 = jlo
  end if

  if ( ncol < jhi ) then
    jhi2 = ncol
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Warning!'
    write ( *, * ) '  Input value of JHI was ',jhi
    write ( *, * ) '  Reset to JHI2 = ',jhi2
  else
    jhi2 = jhi
  end if

  if ( jhi2 < jlo2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PrMatFL - Warning:'
    write ( *, * ) '  Input JLO =      ',jlo, ' JHI= ',jhi
    write ( *, * ) '  Effective JLO2 = ',jlo2,' JHI2=',jhi2
    return
  end if
!
!  Print the individual entries.
!
  write ( *, * ) ' '
  write ( *, * ) 'PrMatFL: Matrix entries'
  write ( *, * ) '  for nodes ILO2 = ',ilo2,' to IHI2=',ihi2
  write ( *, * ) '  and columns JLO2 = ',jlo2,' to JHI2=',jhi2
  write ( *, * ) ' '
!
!  Print columns in groups, going from LLO to LHI.
!
  ncols = 5

  do llo = jlo2,jhi2,ncols

    lhi = min(llo+ncols-1,jhi2)
    write ( *, * ) ' '
    write(*,'(''    Eqn Node'',6(i6,6x))')(l,l = llo,lhi)
    write ( *, * ) ' '
!
!  Compute the index of the equation just before the first
!  equation to be printed.
!
    i = indx(1,ilo2)-1

    do j = ilo2,ihi2

      do k = 1, 3

        if ( 0 < indx(k,j) ) then

          i = i+1

          write(*,'(a2,2i4,6g12.4)')eqn(i),i,j,(a(i,l),l = llo,lhi)

        end if

      end do

    end do

  end do

  return
end
subroutine prpar(iopt,npar,nparb,nparf,par)

!*****************************************************************************80
!
!! PRPAR prints out the current parameters.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!         The number of parameters associated with the position and
!         shape of the bump.
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1. 
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!
!         PAR is the current estimate for the parameters.
!
  implicit none
!
  integer ( kind = 4 ) npar
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iopt(npar)
  character ( len = 6 ) label1
  character ( len = 5 ) label2
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  real ( kind = 8 ) par(npar)
!
  write ( *, * ) ' '

  do i = 1,npar

    if ( iopt(i) == 0) then
      label2 = 'Fixed'
    else
      label2 = 'Free '
    end if

    if ( i <= nparf) then
      label1 = 'Inflow  '
    else if ( nparf+1 <= i.and.i <= nparf+nparb) then
      label1 = 'Bump  '
    else if ( i == nparf+nparb+1) then
      label1 = 'Reynld'
    end if

    write(*,'(i2,2x,a6,2x,a5,2x,g14.6)')i,label1,label2,par(i)

  end do

  return
end
subroutine prsenn(maxcofrb,maxnfl,ncofrb,neqnfl,senfl)

!*****************************************************************************80
!
!! PRSENN prints out the norms of the sensitivities.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  MAXNFL Input, integer ( kind = 4 ) MAXNFL.
!         MAXNFL is the maximum number of equations or coefficients allowed
!         for the full system.  MAXNFL must be used instead of NEQNFL as
!         the leading dimension of certain multi-dimensional arrays.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL.
!         NEQNFL is the number of equations (and coefficients) in the full
!         finite element system.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB.
!         NCOFRB is the number of basis functions used for the
!         reduced basis method.  (The first basis vector is labeled
!         "0").  In this program, that amounts to the number of columns
!         in the matrix RB.  NCOFRB is also the number of reduced basis
!         state equations, and reduced basis coefficients GRB.
!
!  SENFL  Input, real ( kind = 8 ) SENFL(MAXNFL,NCOFRB).
!         SENFL contains the first several order sensitivities of the
!         full solution with respect to the REYNLD parameter.
!
!         SENFL(I,J) contains the J-th sensitivity of the I-th full unknown
!         with respect to REYNLD.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) neqnfl
!
  real ( kind = 8 ) anrmg
  real ( kind = 8 ) dnrm2
  real ( kind = 8 ) enrmg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) itemp
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
!
!  Print the norms of the columns of SENFL.
!
  write ( *, * ) ' '
  write ( *, * ) 'Order    MxNorm     Index    l2 Norm'
  write ( *, * ) ' '

  do i = 1,ncofrb
    itemp = idamax(neqnfl,senfl(1,i),1)
    anrmg = abs(senfl(itemp,i))
    enrmg = dnrm2(neqnfl,senfl(itemp,i),1)
    write(*,'(i6,g14.6,i6,g14.6)')i,anrmg,itemp,enrmg
  end do

  return
end
subroutine prvecfl(eqn,ihi,ilo,indx,neqnfl,np,vec)

!*****************************************************************************80
!
!! PRVECFL prints out some entries of a vector indexed by node number.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IHI,
!  ILO    Input, integer ( kind = 4 ) IHI, ILO.
!         ILO is the first node, and IHI the last node, for which the
!         data should be printed.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K).
! 
!         If INDX(I,J) is positive, then that means that a degree of
!         freedom for variable I (U, V or P) is associated with node
!         J, and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL.
!         NEQNFL is the number of equations in the full system.
!
!  NP     Input, integer ( kind = 4 ) NP.
!         NP is the number of nodes.
!
!  VEC    Input, real ( kind = 8 ) VEC(NEQNFL).
!         VEC contains the vector information to be printed.
!
  implicit none
!
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  character ( len = 2 ) eqn(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi2
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo2
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind = 8 ) vec(neqnfl)
!
  ilo2 = max(1,ilo)
  ihi2 = min(np,ihi)

  if ( ihi2 < ilo2 ) then
    write ( *, * ) ' '
    write ( *, * ) 'PrVecFL - Warning:'
    write ( *, * ) '  Input ILO =      ',ilo, ' IHI= ',ihi
    write ( *, * ) '  Effective ILO2 = ',ilo2,' IHI2=',ihi2
    return
  end if
!
!  Print the individual entries.
!
  write ( *, * ) ' '
  write ( *, * ) 'PrVecFL: Vector entries'
  write ( *, * ) '  for nodes ILO2 = ',ilo2,' to IHI2=',ihi2
  write ( *, * ) ' '
  write ( *, * ) '   Eqn Node      Value'
  write ( *, * ) ' '
!
!  Compute the index of the equation just before the first
!  equation to be printed.
!
  i = indx(1,ilo2)-1

  do j = ilo2,ihi2

    do k = 1, 3

      if ( 0 < indx(k,j) ) then

        i = i+1

        write(*,'(a2,2i5,5g12.4)')eqn(i),i,j,vec(i)

      end if

    end do

  end do

  return
end
subroutine prvecrb ( ihi, ilo, nhi, nlo, vec )
!
!*****************************************************************************80
!
!! PRVECRB prints out entries ILO through IHI of a vector.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IHI,
!  ILO    Input, integer ( kind = 4 ) IHI, ILO.
!         ILO is the lowest, and IHI the highest index to print.
!
!  NHI,
!  NHI    Input, integer ( kind = 4 ) NHI, NLO.
!         NLO is the index of the lowest, and NHI the index of
!         the highest entry of VEC.
!
!  VEC    Input, real ( kind = 8 ) VEC(NLO:NHI).
!         VEC is the vector whose entries are to be printed.
!
  implicit none
!
  integer ( kind = 4 ) nhi
  integer ( kind = 4 ) nlo
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ihi2
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ilo2
  real ( kind = 8 ) vec(nlo:nhi)
!
  if ( nhi < nlo ) then
    write ( *, * ) ' '
    write ( *, * ) 'PrVecRB - Fatal error!'
    write ( *, * ) '  NHI < NLO!'
    write ( *, * ) '  NLO = ',nlo
    write ( *, * ) '  NHI = ',nhi
    stop
  end if

  if ( ilo < nlo) then
    ilo2 = nlo
  else
    ilo2 = ilo
  end if

  if ( nhi < ihi ) then
    ihi2 = nhi
  else
    ihi2 = ihi
  end if

  write ( *, * ) 'Vector entries ILO2 = ',ilo2,' to IHI2=',ihi2
  write ( *, * ) ' '
  do i = ilo2,ihi2
    write(*,'(i6,g14.6)')i,vec(i)
  end do

  return
end
subroutine prxy(ihi,ilo,np,ny,xc,yc)
!
!*****************************************************************************80
!
!! PRXY prints the X and Y coordinates of each node.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IHI,
!  ILO    Input, integer ( kind = 4 ) IHI, ILO, the indices of the last and
!         first nodes to print.
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  NY     Input, integer ( kind = 4 ) NY, the number of nodes in each column.
!
!  XC     Input, real XC(NP), the X coordinates of the nodes.
!
!  YC     Input, real YC(NP), the Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) np
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) ny
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)
!
  write ( *, * ) ' '
  write ( *, * ) 'PrXY:'
  write ( *, * ) '  Print X and Y coordinates of nodes'
  write ( *, * ) ilo,' through ',ihi
  write ( *, * ) ' '
  write ( *, * ) '  Node  Row  Column     X      Y'
  write ( *, * ) ' '

  do i = max(ilo,1),min(ihi,np)
    icol = (i-1)/(2*ny-1) + 1
    irow = i-(icol-1)*(2*ny-1)
    if ( irow == 1) then
      write ( *, * ) ' '
    end if
    write(*,'(3i5,2f12.5)')i,irow,icol,xc(i),yc(i)
  end do

  return
end
subroutine qbf(ielem,in,w,dwdx,dwdy,nelem,node,np,xc,xq,yc,yq)
!
!*****************************************************************************80
!
!! QBF evaluates a particular quadratic basis function at a point
!  in a nonisoparametric element.
!
!      ^
!      |        2
!      |       /|
!   Y  |      4 5
!      |     /  |
!      |    1-6-3
!      |
!      +------------>
!             X
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IELEM  Input, integer ( kind = 4 ) IELEM, the number of the element we are
!         examining.  This will be a value between 1 and NELEM.
!
!  IN     Input, integer ( kind = 4 ) IN, the number of the basis function we
!         want.  This will be a value between 1 and 6.  Functions
!         1 through 3 are associated with corners, 4 though 6
!         with sides.
!
!  W,
!  DWDX,
!  DWDY   Output, real ( kind = 8 ) W, DWDX, DWDY, the value of the
!         IN-th basis  function and its X and Y derivatives, at the
!         given point.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM), contains the numbers
!         of the nodes that make up each element.  Element number
!         I is associated with nodes NODE(1,I) through NODE(6,I).
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  XC     Input, real ( kind = 8 ) XC(NP), the X coordinates of the
!         nodes.
!
!  XQ     Input, real ( kind = 8 ) XQ, the X coordinate of the point
!         where the basis function is to be evaluated.
!
!  YC     Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes
!
!  YQ     Input, real ( kind = 8 ) YQ, the Y coordinate of the point wher
!         the basis function is to be evaluated.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) c
  real ( kind = 8 ) d
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) in
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) in3
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) s
  real ( kind = 8 ) t
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq
!
!  Case 1: We are inquiring about a basis function associated
!  with a corner.
!
!  Notice that the basis function W is zero exactly if
!  T is 0 or T is 1/2.
!
!  IN1, IN2, and IN3 are the local node numbers of the three
!  corner nodes, and I1, I2 and I3 are the corresponding
!  global node numbers, which are used to look up the X and
!  Y coordinates of the nodes.
!
  if ( 1 <= in.and.in <= 3) then

    in1 = in
    in2 = mod(in,3)+1
    in3 = mod(in+1,3)+1

    i1 = node(in1,ielem)
    i2 = node(in2,ielem)
    i3 = node(in3,ielem)

    d = (xc(i2)-xc(i1))*(yc(i3)-yc(i1)) -(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

    if ( d == 0.0D+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'QBF - Fatal error!'
      write ( *, * ) '  D = 0'
      write ( *, * ) '  Element IELEM = ',ielem
      write ( *, * ) '  I1, XC(I1), YC(I1) = ',i1,xc(i1),yc(i1)
      write ( *, * ) '  I2, XC(I2), YC(I2) = ',i2,xc(i2),yc(i2)
      write ( *, * ) '  I3, XC(I3), YC(I3) = ',i3,xc(i3),yc(i3)
      stop
    end if

    t = 1.0D+00 + ( (xq -xc(i1))*(yc(i2)-yc(i3))+(xc(i3)-xc(i2))*(yq -yc(i1)) )/d

    w = t*(2.0D+00 *t-1.0D+00 )

    dwdx = (yc(i2)-yc(i3))*(4.0D+00 *t-1.0D+00 )/d
    dwdy = (xc(i3)-xc(i2))*(4.0D+00 *t-1.0D+00 )/d
!
!  Case 2: We are inquiring about a basis function associated
!  with a midpoint.
!
  else if ( 4 <= in .and. in <= 6 ) then

    in1 = in-3
    in2 = mod(in-3,3)+1
    in3 = mod(in-2,3)+1

    i1 = node(in1,ielem)
    i2 = node(in2,ielem)
    i3 = node(in3,ielem)

    d =     (xc(i2)-xc(i1))*(yc(i3)-yc(i1))-(xc(i3)-xc(i1))*(yc(i2)-yc(i1))

    if ( d == 0.0D+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'QBF - Fatal error!'
      write ( *, * ) '  D = 0'
      write ( *, * ) '  Element IELEM = ',ielem
      write ( *, * ) '  I1, XC(I1), YC(I1) = ',i1,xc(i1),yc(i1)
      write ( *, * ) '  I2, XC(I2), YC(I2) = ',i2,xc(i2),yc(i2)
      write ( *, * ) '  I3, XC(I3), YC(I3) = ',i3,xc(i3),yc(i3)
      stop
    end if

    c =     (xc(i3)-xc(i2))*(yc(i1)-yc(i2))-(xc(i1)-xc(i2))*(yc(i3)-yc(i2))

    if ( c == 0.0D+00 ) then
      write ( *, * ) ' '
      write ( *, * ) 'QBF - Fatal error!'
      write ( *, * ) '  C = 0'
      write ( *, * ) '  Element IELEM = ',ielem
      write ( *, * ) '  I1, XC(I1), YC(I1) = ',i1,xc(i1),yc(i1)
      write ( *, * ) '  I2, XC(I2), YC(I2) = ',i2,xc(i2),yc(i2)
      write ( *, * ) '  I3, XC(I3), YC(I3) = ',i3,xc(i3),yc(i3)
      stop
    end if

    t = 1.0D+00 + ( (xq-xc(i1))*(yc(i2)-yc(i3))+(xc(i3)-xc(i2))*(yq    -yc(i1)) )/d

    s = 1.0D+00 + ( (xq-xc(i2))*(yc(i3)-yc(i1))+(xc(i1)-xc(i3))*(yq    -yc(i2)) )/c

    w = 4.0D+00 * s*t
    dwdx = 4.0D+00 * ((yc(i3)-yc(i1))*t/c + (yc(i2)-yc(i3))*s/d)
    dwdy = 4.0D+00 * ((xc(i1)-xc(i3))*t/c + (xc(i3)-xc(i2))*s/d)

  else

    write ( *, * ) ' '
    write ( *, * ) 'QBF - Fatal error!'
    write ( *, * ) '  Request for basis function IN = ',in
    write ( *, * ) '  but IN must be between 1 and 6.'
    stop

  end if

  return
end
subroutine refbsp(q,dqdx,dqdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)
!
!*****************************************************************************80
!
!! REFBSP evaluates one of the three linear basis functions,
!  and its X and Y derivatives, at a particular point (X,Y)
!  in a particular element, by referring to the corresponding
!  points (XSI,ETA) in the reference triangle.
!
!  It is assumed that we already know the value of the jacobian
!  of the isoparametric transformation between the (XSI, ETA) and
!  (X, Y) spaces.  The four entries of the jacobian are
!  symbolically named DETADX, DETADY, DXSIDX and DXSIDY, and
!  we know that the jacobian gives us the following relation
!  between derivatives with respect to XSI and ETA, and derivatives
!  with respect to X and Y:
!
!    dF/dX = dF/dXsi dXsi/dX + dF/dEta dEta/dX
!    dF/dY = dF/dXsi dXsi/dY + dF/dEta dEta/dY
!
!  Here is a graph of the (XSI, ETA) reference triangle we will
!  use.
!
!        ^
!        |
!      1 +        2
!        |       /|
!  ETA   |      / |
!        |     /  |
!      0 +    1---3
!        |
!        +----+---+--->
!             0   1
!
!              XSI
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  Q,
!  DQDX,
!  DQDY   Output, real ( kind = 8 ) Q, DQDX, DQDY, the value of the basis
!         function, and its derivatives with respect to X and Y, at
!         the point (ETA,XSI).
!
!  DETADX,
!  DETADY Input, real ( kind = 8 ) DETADX, DETADY, the partial derivative
!         d ETA/d X and d ETA/d Y at (ETA,XSI).
!
!  IQ     Input, integer ( kind = 4 ) IQ, the local node number, between 1 and
!         3, whose basis function is being evaluated.
!
!  DXSIDX,
!  DXSIDY Input, real ( kind = 8 ) DXSIDX, DXSIDY, the partial derivative
!         d XSI/d X and d XSI/d Y at (ETA,XSI).
!
!  ETA,
!  XSI    Input, real ( kind = 8 ) ETA, XSI, the local coordinates of the
!         at which the basis information is desired.
!
  implicit none
!
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dqdeta
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdxsi
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  integer ( kind = 4 ) iq
  real ( kind = 8 ) q
  real ( kind = 8 ) xsi
!
!  Refuse to evaluate the basis functions for arguments (XSI,ETA)
!  that lie outside the reference triangle.
!
  if ( 1.0D+00 < xsi ) then
    write ( *, * ) ' '
    write ( *, * ) 'REFBSP - Fatal error!'
    write ( *, * ) '  XSI must be less than or equal to 1.'
    write ( *, * ) '  Input XSI is ',xsi
    stop
  end if

  if ( xsi < eta  ) then
    write ( *, * ) ' '
    write ( *, * ) 'REFBSP - Fatal error!'
    write ( *, * ) '  ETA must be less or equal to XSI.'
    write ( *, * ) '  Input XSI, ETA = ',xsi,eta
    stop
  end if

  if ( eta < 0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'REFBSP - Fatal error!'
    write ( *, * ) '  ETA must be 0 or greater.'
    write ( *, * ) '  Input ETA = ',eta
    stop
  end if
!
  if ( iq == 1) then
    q = 1.0D+00 -xsi
    dqdxsi = -1.0D+00
    dqdeta =  0.0D+00
  else if ( iq == 2) then
    q = eta
    dqdxsi = 0.0D+00
    dqdeta = 1.0D+00
  else if ( iq == 3) then
    q = xsi-eta
    dqdxsi = 1.0D+00
    dqdeta = -1.0D+00
  else if ( 4 <= iq .and. iq <= 6 ) then
    q = 0.0D+00
    dqdxsi = 0.0D+00
    dqdeta = 0.0D+00
  else
    write ( *, * ) 'RefBSP - Fatal error!'
    write ( *, * ) '  Request for basis function IQ = ',iq
    write ( *, * ) '  but IQ must be between 1 and 6.'
    stop
  end if

  dqdx = dqdxsi*dxsidx+dqdeta*detadx
  dqdy = dqdxsi*dxsidy+dqdeta*detady

  return
end
subroutine refqbf(w,dwdx,dwdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)
!
!*****************************************************************************80
!
!! REFQBF evaluates one of the six quadratic basis functions,
!  and its X and Y derivatives, at a particular point in a
!  particular element, by referring to the reference triangle.
!
!  The point we are interested in is referred to by its coordinates
!  in the reference triangle.  That is, we are given coordinates
!  (XSI, ETA), even though, physically, we are interested
!  in points in (X, Y) space.
!
!  It is assumed that we already know the value of the jacobian
!  of the isoparametric transformation between the (XSI, ETA) and
!  (X, Y) spaces.  The four entries of the jacobian are
!  symbolically named DETADX, DETADY, DXSIDX and DXSIDY, and
!  we know that the jacobian gives us the following relation
!  between derivatives with respect to XSI and ETA, and derivatives
!  with respect to X and Y:
!
!    d F(X,Y)/dX     (d XSI/dX  d ETA/dX )   ( d F(XSI, ETA)/d XSI )
!    d F(X,Y)/dY  =  (d XSI/dY  d ETA/dY ) * ( d F(XSI, ETA)/d ETA )
!
!  Here is a graph of the (XSI, ETA) reference triangle we will
!  use.
!
!        ^
!        |
!      1 +        2
!        |       /|
!  ETA   |      4 5
!        |     /  |
!      0 +    1-6-3
!        |
!        +----+---+--->
!             0   1
!
!              XSI
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  W,
!  DWDX,
!  DWDY   Output, real ( kind = 8 ) W, DWDX, DWDY, the value of the basis
!         function, and its derivatives with respect to X and Y, at
!         the point (XSI,ETA).
!
!  DETADX,
!  DETADY Input, real ( kind = 8 ) DETADX, DETADY, the partial derivative
!         d ETA/d X and d ETA/d Y at (XSI,ETA).
!
!  DXSIDX,
!  DXSIDY Input, real ( kind = 8 ) DXSIDX, DXSIDY, the partial derivative
!         d XSI/d X and d XSI/d Y at (XSI,ETA).
!
!  ETA    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!  IQ     Input, integer ( kind = 4 ) IQ, the local node number, between 1 and
!         6, whose basis function is being evaluated.
!
!  XSI    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
  implicit none
!
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dwdeta
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdxsi
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  integer ( kind = 4 ) iq
  real ( kind = 8 ) w
  real ( kind = 8 ) xsi
!
!  Refuse to evaluate the basis functions for arguments (XSI,ETA)
!  that lie outside the reference triangle.
!
  if ( 1.0D+00 < xsi ) then
    write ( *, * ) ' '
    write ( *, * ) 'REFQBF - Fatal error!'
    write ( *, * ) '  XSI must be less than or equal to 1.'
    write ( *, * ) '  Input XSI is ',xsi
    stop
  end if

  if ( xsi < eta ) then
    write ( *, * ) ' '
    write ( *, * ) 'REFQBF - Fatal error!'
    write ( *, * ) '  ETA must be less or equal to XSI.'
    write ( *, * ) '  Input XSI, ETA = ',xsi,eta
    stop
  end if

  if ( eta < 0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'REFQBF - Fatal error!'
    write ( *, * ) '  ETA must be 0 or greater.'
    write ( *, * ) '  Input ETA = ',eta
    stop
  end if
!
!  Evaluate W, the quadratic basis function.
!  Evaluate DWDXSI and DWDETA, the partial derivatives d W/d XSI
!  and d W/d ETA.
!
!  Basis 1 is zero if XSI = 0.5 or XSI=1.
!
  if ( iq == 1) then
    w =  (2.0D+00 *xsi-1.0D+00 ) * (xsi-1.0)
    dwdxsi = -3.0D+00 + 4.0D+00 *xsi
    dwdeta = 0.0D+00
!
!  Basis 2 is zero if ETA = 0 or ETA=0.5.
!
  else if ( iq == 2) then
    w =  eta * (2.0D+00 *eta-1.0D+00 )
    dwdxsi = 0.0D+00
    dwdeta = -1.0D+00 + 4.0D+00 *eta
!
!  Basis 3 is zero if XSI = ETA, or XSI=ETA+0.5
!
  else if ( iq == 3) then
    w =  (xsi-eta) * (2.0D+00 *xsi-2.0D+00 *eta-1.0D+00 )
    dwdxsi = -1.0D+00 + 4.0D+00 *xsi-4.0D+00 *eta
    dwdeta = 1.0D+00 -4.0D+00 *xsi+4.0D+00 *eta
!
!  Basis 4 is zero if ETA = 0 or XSI=1.
!
  else if ( iq == 4) then
    w =  4.0D+00 * eta * (1.0D+00 -xsi)
    dwdxsi = -4.0D+00 *eta
    dwdeta = 4.0D+00 -4.0D+00 *xsi
!
!  Basis 5 is zero if ETA = 0 or XSI=ETA.
!
  else if ( iq == 5) then
    w = 4.0D+00 * eta * (xsi-eta)
    dwdxsi = 4.0D+00 *eta
    dwdeta = 4.0D+00 *xsi-8.0D+00 *eta
!
!  Basis 6 is zero if XSI = ETA or XSI=1.
!
  else if ( iq == 6) then
    w = 4.0D+00 * (xsi-eta) * (1.0D+00 -xsi)
    dwdxsi = 4.0D+00 -8.0D+00 *xsi+4.0D+00 *eta
    dwdeta = -4.0D+00 + 4.0D+00 * xsi
!
!  Stop if we were given an unexpected value of IQ.
!
  else
    write ( *, * ) ' '
    write ( *, * ) 'RefQBF - Fatal error!'
    write ( *, * ) '  A basis function index must be between 1 and 6,'
    write ( *, * ) '  but you input the value IQ = ',iq
    stop
  end if
!
!  Convert the d W/d XSI and d W/d ETA derivatives to d W/d X
!  and d W/d Y.
!
  dwdx = dwdxsi*dxsidx + dwdeta*detadx
  dwdy = dwdxsi*dxsidy + dwdeta*detady

  return
end
function s_eqi ( s1, s2 )
!
!*****************************************************************************80
!
!! S_EQI is a case insensitive comparison of two strings for equality.
!
!
!  Example:
!
!    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S1, S2, the strings to compare.
!
!    Output, logical S_EQI, the result of the comparison.
!
  character c1
  character c2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) len1
  integer ( kind = 4 ) len2
  integer ( kind = 4 ) lenc
  logical s_eqi
  character ( len = * ) s1
  character ( len = * ) s2
!
  len1 = len ( s1 )
  len2 = len ( s2 )
  lenc = min ( len1, len2 )

  s_eqi = .false.

  do i = 1, lenc

    c1 = s1(i:i)
    c2 = s2(i:i)
    call ch_cap ( c1 )
    call ch_cap ( c2 )

    if ( c1 /= c2 ) then
      return
    end if

  end do

  do i = lenc + 1, len1
    if ( s1(i:i) /= ' ' ) then
      return
    end if
  end do

  do i = lenc + 1, len2
    if ( s2(i:i) /= ' ' ) then
      return
    end if
  end do

  s_eqi = .true.

  return
end
subroutine setban ( indx, ldafl, nelem, nlband, node, np )
!
!*****************************************************************************80
!
!! SETBAN computes NLBAND, the lower band width of the Jacobian matrix
!  stored in LINPACK general band storage format.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  LDAFL  Input, integer ( kind = 4 ) LDAFL, the first dimension of the matrix AFL.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NLBAND Output, integer ( kind = 4 ) NLBAND.
!
!         The lower bandwidth of the matrix A.  The zero structure of A
!         is assumed to be symmetric, and so NLBAND is also the upper
!         bandwidth of A. 
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM).
!
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipp
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iqq
  integer ( kind = 4 ) iuk
  integer ( kind = 4 ) iukk
  integer ( kind = 4 ) j
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,nelem)
!
  nlband = 0

  do ielem = 1,nelem
    do iq = 1,6
      ip = node(iq,ielem)
      do iuk = 1,3
        i = indx(iuk,ip)
        if ( 0 < i ) then
          do iqq = 1, 6
            ipp = node(iqq,ielem)
            do iukk = 1,3
              j = indx(iukk,ipp)
              if ( 0 < j ) then
                if ( nlband < j-i ) then
                  nlband = j-i
                end if
              end if
            end do
          end do
        end if
      end do
    end do
  end do

  if ( ldafl < 3*nlband+1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SetBan - Fatal error!'
    write ( *, * ) '  Not enough room for matrix!'
    write ( *, * ) '  Number of rows needed is ',3*nlband+1
    write ( *, * ) '  The maximum allowed is LDAFL =  ',ldafl
    stop
  end if

  return
end
subroutine setgeo(area,etaq,gridx,gridy,ibs,isotri,nelem,node,nodelm,np, &
  npar,nparb,nparf,nx,ny,par,phifl,region,splbmp,taubmp,wquad,xbl,xbr,xc, &
  xquad,xrange,xsiq,ybl,ybr,yc,yquad,yrange)
!
!*****************************************************************************80
!
!! SETGEO is given a set of flow parameters in PAR, and an
!  approximate solution vector G, and proceeds to set up the
!  constraints associated with PAR, and use Newton iteration
!  to correct G to a solution that satisfies the constraints
!  to within some tolerance.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) AREA(3,NELEM).
!    AREA contains a common factor multiplying the term associated
!    with a quadrature point in a given element, namely,
!      AREA(IQUAD,IELEM) = Ar(IELEM) * WQUAD(IQUAD)
!    or, if the element is isoperimetric,
!      AREA(IQUAD,IELEM) = DET * Ar(IELEM) * WQUAD(IQUAD)
!    Here Ar(IELEM) represents the area of element IELEM.
!
!    Input, real ( kind = 8 ) ETAQ(3).
!    The "Eta" coordinates of the quadrature points.
!
!    Input, integer ( kind = 4 ) IBS.
!    1, the bump is modeled by C0 linear splines.
!    2, the bump is modeled by C0 quadratic splines.
!
!    Input, integer ( kind = 4 ) ISOTRI(NELEM).
!    0, the element is NOT isoparametric, and the nodes never move.
!    That means that the quadrature points are only computed once.
!    1, the element is NOT isoparametric, but the nodes may move.
!    Quadrature point locations must be updated on each step.
!    This could occur for elements above, but not touching, the bump.
!    2, the element is isoparametric.
!
!    Input, integer ( kind = 4 ) NELEM, the number of elements.
!  
!    Input, integer ( kind = 4 ) NODE(6,NELEM).
!    NODE(I,J) contains, for an element J, the global node index of
!    the element node whose local number is I.
!    The local ordering of the nodes is suggested by this diagram:
!
!          2
!         /|
!        4 5
!       /  |
!      1-6-3
!
!    Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!    element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!    Input, integer ( kind = 4 ) NPAR.
!    The number of parameters.  NPAR = NPARF + NPARB + 1.
!    The parameters control the shape of the inflow,
!    the shape of the bump obstacle, and the strength of the
!    flow.
!
!    Input, integer ( kind = 4 ) NPARB.
!    The number of parameters associated with the position and
!    shape of the bump.
!    Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!    Input, integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1. 
!
!    Input, integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!    Roughly speaking, NX (or 2*NX) is the number of elements along
!    a line in the X direction.
!
!    Input, integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!    Roughly speaking, NY (or 2*NY) is the number of elements along
!    a line in the Y direction.
!
!    Input, real ( kind = 8 ) PAR(NPAR).
!    PAR is the current estimate for the parameters.
!
!    Output, real ( kind = 8 ) PHIFL(3,6,10,NELEM). 
!    PHIFL contains the value of a finite element basis function, its
!    derivative, or other information, evaluated at the quadrature
!    points.
!    The meaning of the entry PHIFL(I,J,K,L) is as follows. 
!    For the quadrature point I, and basis function J, in element L,
!    PHIFL(I,J,K,L) represents the value of:
!      K =  1, W, the finite element basis function for velocities;
!      K =  2, dWdX, the X derivative of W;
!      K =  3, dWdY, the Y derivative of W;
!      K =  4, Q, the finite element basis function for pressures;
!      K =  5, dQdX, the X derivative of Q;
!      K =  6, dQdY, the Y derivative of Q;
!      K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!      K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!      K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!      K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!    In particular, PHIFL(I,J,K,L) is the value of the quadratic
!    basis function W associated with local node J in element L,
!    evaluated at quadrature point I.
!    Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!    since there are only three linear basis functions.
!
!    Input, character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottom with tangential velocity specifications there.
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!    Output, real ( kind = 8 ) SPLBMP(NPARB+2).
!    SPLBMP contains the spline coefficients for the bump.
!
!    Output, real ( kind = 8 ) SPLFLO(NPARF+2).
!    SPLFLO contains the spline coefficients for the inflow.
!
!    Output, real ( kind = 8 ) TAUBMP(NPARB+2).
!    TAUBMP contains the location of the spline abscissas for
!    the bump.  There are NPARB+2 of them, because the end values
!    of the spline are constrained to have particular values.
!
!    Output, real ( kind = 8 ) TAUFLO(NPARF+2).
!    TAUFLO contains the location of the spline abscissas for
!    the inflow.  There are NPARF+2 of them, because the end
!    values of the spline are constrained to have particular
!    values.
!
!    Input, real ( kind = 8 ) WQUAD(3), the weights for Gaussian
!    quadrature.
!
!    Input, real ( kind = 8 ) XBL, the X coordinate of the left corner
!    of the bump.
!
!    Input, real ( kind = 8 ) XBR, the X coordinate of the right corner
!    of the bump.
!
!    Input, real ( kind = 8 ) XC(NP).
!    The X coordinates of the nodes.
!
!    Input, real ( kind = 8 ) XQUAD(3,NELEM).
!    The X coordinates of the quadrature points for each element.
!
!    Input, real ( kind = 8 ) XSIQ(3).
!    The "Xsi" coordinates of the quadrature points.
!
!    Input, real ( kind = 8 ) YBL, the Y coordinate of the left corner
!    of the bump.
!
!    Input, real ( kind = 8 ) YBR, the Y coordinate of the right corner
!    of the bump.
!
!    Input, real ( kind = 8 ) YC(NP).
!    The Y coordinates of the nodes.
!
!    Input, real ( kind = 8 ) YQUAD(3,NELEM).
!    The Y coordinates of the quadrature points for each element.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
!
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) etaq(3)
  character ( len = 20 ) gridx
  character ( len = 20 ) gridy
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nodelm(np)
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  character ( len = 20 ) region
  real ( kind = 8 ) splbmp(nparb+2)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(3,nelem)
  real ( kind = 8 ) xrange
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(3,nelem)
  real ( kind = 8 ) yrange
!
!  Set the spline coefficients for the bump.
!
  call bmpspl(npar,nparb,nparf,par,splbmp,taubmp,xbl,xbr,ybl,ybr)
!
!  Set the X and Y coordinates of the nodes that form the grid.
!
  call setxy(gridx,gridy,ibs,np,nparb,nx,ny,region,splbmp,taubmp, &
    xbl,xbr,xc,xrange,ybl,ybr,yc,yrange)
!
!  Set the quadrature points, which move every step if there
!  are bump parameters.
!
  call setq3(area,etaq,isotri,nelem,node,np,wquad,xc,xquad,xsiq,yc,yquad)
!
!  Set the value of the basis functions at all quadrature points.
!
  call setpfl(area,etaq,isotri,nelem,node,np,phifl,xc,xquad,xsiq,yc,yquad)
!
!  Set the NODELM array.
!
  do ip = 1,nelem
    nodelm(ip) = 0
  end do

  do ielem = 1,nelem
    do j = 1,6
      ip = node(j,ielem)
      nodelm(ip) = ielem
    end do
  end do

  return
end
subroutine setlog(eqn,hx,hy,ibump,indx,isotri,ldafl,maxelm, &
  maxnfl,maxnp,nelem,neqnfl,nlband,node,np,nprof,nx,ny, &
  region,xbl,xbr,xprof,xrange,ybr,yrange)
!
!*****************************************************************************80
!
!! SETLOG determines some data that depends on the user input.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  EQN    Output, character ( len = 2 ) EQN(NEQNFL).
!         EQN records the "type" of each equation that will be generated, and
!         which is associated with an unknown.  Note that most boundary
!         conditions do not result in an equation.  The current values are:
!
!         'U'  The horizontal momentum equation.
!         'UB' The condition U = 0 applied at a node on the bump.
!         'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!         'UW' The condition U = 0 applied at a node on a fixed wall.
!
!         'V'  The vertical momentum equation.
!         'VB' The condition V = 0 applied at a node on the bump.
!         'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!         'VW' The condition V = 0 applied at a node on a fixed wall.
!
!         'P'  The continuity equation.
!         'PB' The condition P = 0 applied at (XMAX,YMAX).
!
!  HX     Output, real ( kind = 8 ) HX.
!         HX is the nominal spacing between nodes in the X direction.
!
!  HY     Output, real ( kind = 8 ) HY.
!         HY is the nominal spacing between nodes in the Y direction.
!
!  NELEM  Output, integer ( kind = 4 ) NELEM.
!         NELEM is the number of elements.
!         NELEM can be determined as 2*(NX-1)*(NY-1).
!
!  NP     Output, integer ( kind = 4 ) NP.
!         NP is the number of nodes used to define the finite element mesh.
!         Typically, the mesh is generated as a rectangular array, with
!         an odd number of nodes in the horizontal and vertical directions.
!         The formula for NP is NP = (2*NX-1)*(2*NY-1).
!
!  NX     Input, integer ( kind = 4 ) NX.
!
!         NX controls the spacing of nodes and elements in
!         the X direction.  There are 2*NX-1 nodes along various
!         lines in the X direction.
!
!         Roughly speaking, NX (or 2*NX) is the number of elements along
!         a line in the X direction.
!
!  NY     Input, integer ( kind = 4 ) NY.
!
!         NY controls the spacing of nodes and elements in
!         the Y direction.  There are 2*NY-1 nodes along various
!         lines in the Y direction.
!
!         Roughly speaking, NY (or 2*NY) is the number of elements along
!         a line in the Y direction.
!
!  XPROF  Output, real ( kind = 8 ) XPROF.
!
!         The X coordinate at which the profile is measured.  This
!         value should be a grid value!
!
!  XRANGE Input, real ( kind = 8 ) XRANGE.
!         The total width of the region.
!
!  YRANGE Input, real ( kind = 8 ) YRANGE.
!         YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) maxnp
  integer ( kind = 4 ) ny
!
  character ( len = 2 ) eqn(maxnfl)
  real ( kind = 8 ) hx
  real ( kind = 8 ) hy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) indx(3,maxnp)
  integer ( kind = 4 ) isotri(maxelm)
  integer ( kind = 4 ) itemp
  integer ( kind = 4 ) ldafl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nlband
  integer ( kind = 4 ) node(6,maxelm)
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) nx
  character ( len = 20 ) region
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yrange
!
  nelem = 2*(nx-1)*(ny-1)

  write ( *, * ) ' '
  write ( *, * ) 'SetLog - Note:'
  write ( *, * ) '  Number of elements, NELEM = ',nelem

  np = (2*nx-1)*(2*ny-1)

  write ( *, * ) '  Number of nodes, NP = ',np

  if ( xprof < 0.0 .or. xrange < xprof ) then
    write ( *, * ) ' '
    write ( *, * ) 'SetLog - Fatal error!'
    write ( *, * ) '  XPROF lies outside of XRANGE.'
    write ( *, * ) '  XPROF = ',xprof
    write ( *, * ) '  XRANGE = ',xrange
    stop
  end if

  if ( 1 < nx ) then
    hx = xrange/(2.0D+00 *dble(nx-1))
  else
    write ( *, * ) ' '
    write ( *, * ) 'SetLog - Fatal error!'
    write ( *, * ) '  NX = ',nx
    stop
  end if

  if ( 1 < ny ) then
    hy = yrange/(2.0D+00 *dble(ny-1))
  else
    write ( *, * ) ' '
    write ( *, * ) 'SetLog - Fatal error!'
    write ( *, * ) '  NY = ',ny
    stop
  end if

  write ( *, * ) '  X nodal spacing is HX = ',hx
  write ( *, * ) '  Y nodal spacing is HY = ',hy

!
!  Set the logical NODE array.
!
  call setnod(eqn,ibump,indx,isotri,maxnfl,nelem,neqnfl,node, &
    np,nx,ny,region,xbl,xbr,xrange,ybr,yrange)

  write ( *, * ) '  The number of unknowns is NEQNFL = ',neqnfl

!
!  Set the location of the profile nodes.
!
  itemp = nint((2.0D+00 *dble(nx-1)*xprof)/xrange)
  do i = 1,2*ny-1
    nprof(i) = itemp*(2*ny-1)+i
  end do

  write ( *, * ) ' '
  write ( *, * ) '  Profile nodes extend from ',nprof(1)
  write ( *, * ) '  to ',nprof(2*ny-1)
!
!  Get the matrix bandwidth.
!
  write ( *, * ) ' '
  write ( *, * ) '  Maximum full matrix rows LDAFL = ',ldafl

  call setban(indx,ldafl,nelem,nlband,node,np)

  write ( *, * ) '  Lower bandwidth NLBAND =     ',nlband
  write ( *, * ) '  Required matrix rows 3*NLBAND+1 =  ',3*nlband+1

  return
end
subroutine setnod(eqn,ibump,indx,isotri,maxnfl,nelem,neqnfl,node, &
  np,nx,ny,region,xbl,xbr,xrange,ybr,yrange)
!
!*****************************************************************************80
!
!! SETNOD assigns numbers to the nodes and elements, decides which
!  elements shall be isoparametric, (ISOTRI) and assigns six nodes
!  to each (NODE). 
!
!
!  Discussion:
!
!    The routine associates global unknown indices with each node (INDX), and
!    computes the total number of unknowns and equations (NEQNFL), and
!    compares that to the maximum allowed value, MAXNFL.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  EQN   
!    Output, character ( len = 2 ) EQN(NEQNFL).
!    EQN records the "type" of each equation that will be generated, and
!    which is associated with an unknown.  Note that most boundary
!    conditions do not result in an equation.  The current values are:
!
!    'U'  The horizontal momentum equation.
!    'UB' The condition U = 0 applied at a node on the bump.
!    'UI' The condition U = UInflow(Y,Lambda) at the inflow.
!    'UW' The condition U = 0 applied at a node on a fixed wall.
!
!    'V'  The vertical momentum equation.
!    'VB' The condition V = 0 applied at a node on the bump.
!    'VI' The condition V = VInflow(Y,Lambda) at the inflow.
!    'VW' The condition V = 0 applied at a node on a fixed wall.
!
!    'P'  The continuity equation.
!    'PB' The condition P = 0 applied at (XMAX,YMAX).
!
!  IBUMP 
!    Input, integer ( kind = 4 ) IBUMP. 
!
!    IBUMP determines where isoparametric elements will be used.
!
!    0, no isoparametric elements will be used. 
!       Midside nodes of nonisoparametric elements above the
!       bump will be recomputed so that the sides are straight.
!
!    1, isoparametric elements will be used only for the
!       elements which directly impinge on the bump.
!       Midside nodes of nonisoparametric elements above the
!       bump will be recomputed so that the sides are straight.
!
!    2, isoparametric elements will be used for all elements
!       which are above the bump.  All nodes above the bump
!       will be equally spaced in the Y direction.
!
!    3, isoparametric elements will be used for all elements.
!       All nodes above the bump will be equally spaced in
!       the Y direction.
!
!  INDX  
!    Output, integer ( kind = 4 ) INDX(3,NP). 
!
!    INDX(I,J) contains, for each node J, the index of U, V and P at
!    that node, or 0 or a negative value.
!
!    If K = INDX(I,J) is positive, then the value of the degree
!    of freedom is stored in the solution vector entry GFL(K),
!    and an equation will be generated to determine its value.
!
!    If INDX(I,J) is not positive, then no equation is
!    generated to determine for variable I at node J, either because
!    the variable is specified in some other way, or because
!    (in the case of pressure), there is no coefficient associated
!    with that node.
!
!  ISOTRI
!    Output, integer ( kind = 4 ) ISOTRI(NELEM).
!
!    0, the element is NOT isoparametric.
!
!    1, the element is isoparametric.
!
!  MAXNFL
!    Input, integer ( kind = 4 ) MAXNFL.
!    MAXNFL is the maximum number of equations that
!    the program can handle.
!
!  NELEM 
!    Input, integer ( kind = 4 ) NELEM.
!    NELEM is the number of elements.
!
!  NEQNFL
!    Output, integer ( kind = 4 ) NEQNFL.
!    NEQNFL is the number of finite element equations used
!    to define the horizontal and vertical velocities and the
!    pressure.
!
!  NODE  
!    Output, integer ( kind = 4 ) NODE(6,NELEM).
!    NODE contains the numbers of the nodes that make up each element.
!    Element number I is associated with nodes NODE(1,I) through NODE(6,I).
!
!  NP    
!    Input, integer ( kind = 4 ) NP.
!    NP is the number of nodes.
!
!  NX   
!    Input, integer ( kind = 4 ) NX.
!
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!
!    Roughly speaking, NX (or 2*NX) is the number of elements along
!    a line in the X direction.
!
!  NY    
!    Input, integer ( kind = 4 ) NY.
!
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!
!    Roughly speaking, NY (or 2*NY) is the number of elements along
!    a line in the Y direction.
!
!  REGION
!    Input, character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottom, with tangential velocity specifications
!    there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!  XBL   
!    Input, real ( kind = 8 ) XBL.
!    The X coordinate of the left corner of the bump.
!
!  XBR   
!    Input, real ( kind = 8 ) XBR.
!    The X coordinate of the right corner of the bump.
!
!  XRANGE
!    Input, real ( kind = 8 ) XRANGE.
!    The total width of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  character ( len = 2 ) eqn(maxnfl)
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) icol2
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) irow2
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jrow
  logical s_eqi
  integer ( kind = 4 ) nbleft
  integer ( kind = 4 ) nbrite
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  character ( len = 20 ) region
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yrange
!
!  For the channel problem only, compute the global node numbers that will be
!  assigned to the beginning and ending of the bump.  These numbers are
!  only used to determine which elements are isoparametric.
!
  if ( s_eqi ( region,'channel')) then

    nbleft = nint(xbl*(2*nx-2)/xrange)+1
    nbrite = nint(xbr*(2*nx-2)/xrange)+1

  end if
!
!  For the step problem only, determine the horizontal and vertical
!  mesh lines that will pass through the point (XBR,YBR).  Here, we are
!  assuming the step is vertical, and nonnegative!
!
  if ( s_eqi ( region,'step')) then

    jcol = 2*nint(xbr*(nx-1)/xrange)+1
    jrow = 2*nint(ybr*(ny-1)/yrange)+1

  end if
!
!  Consider each of the NP nodes, which logically lie in an MX by MY
!  rectangular array.  A pair of new elements must be generated every
!  time we reach a node that lies in an odd row and column, (except for
!  the top row, and last column, of course).  At every node, we
!  will have to decide how many equations to generate.
!
  ielem = 0
  neqnfl = 0

  do ip = 1,np
!
!  Determine the row and column of this node, and also whether each
!  of these quantities is odd or even.
!
    icol = ((ip-1)/(2*ny-1))+1
    irow = mod((ip-1),2*ny-1)+1

    icol2 = mod(icol,2)
    irow2 = mod(irow,2)
!
!  If both the row and the column are odd, and we're not in the last
!  column or top row, then we can define two new triangular elements
!  based at the node.
!
!  Given the following arrangement of nodes, for instance:
!
!    05 10 15 20 25
!    04 09 14 19 24
!    03 08 13 18 23
!    02 07 12 17 22
!    01 06 11 16 21
!
!  when we arrive at node 13, we will define
!
!    element 7: (25, 13, 15, 19, 14, 20)
!    element 8: (13, 25, 23, 19, 24, 18)
!
    if ( (irow2 == 1.and.icol2.eq.1).and. &
         (icol /= 2*nx-1).and.(irow.ne.2*ny-1)) then

      ielem = ielem+1

      node(1,ielem) = ip+2*(2*ny-1)+2
      node(2,ielem) = ip
      node(3,ielem) = ip+2
      node(4,ielem) = ip+(2*ny-1)+1
      node(5,ielem) = ip+1
      node(6,ielem) = ip+(2*ny-1)+2
!
!  Determine if the elements are isoparametric.
!
      if ( s_eqi ( region,'cavity')) then

        isotri(ielem) = 0

      else if ( s_eqi ( region,'cavity2')) then

        isotri(ielem) = 0

      else if ( s_eqi ( region,'channel')) then

        if ( ibump == 0) then

          if ( nbleft <= icol .and. icol < nbrite) then
            isotri(ielem) = 1
          else
            isotri(ielem) = 0
          end if

        else if ( ibump == 1) then

          if ( nbleft <= icol .and. icol < nbrite ) then
            isotri(ielem) = 1
          else
            isotri(ielem) = 0
          end if

        else if ( ibump == 2) then

          if ( nbleft <= icol .and. icol < nbrite ) then
            isotri(ielem) = 2
          else
            isotri(ielem) = 0
          end if

        else

          isotri(ielem) = 2

        end if

      else if ( s_eqi ( region,'step')) then

        isotri(ielem) = 0

      end if

      ielem = ielem+1

      node(1,ielem) = ip
      node(2,ielem) = ip+2*(2*ny-1)+2
      node(3,ielem) = ip+2*(2*ny-1)
      node(4,ielem) = ip+(2*ny-1)+1
      node(5,ielem) = ip+2*(2*ny-1)+1
      node(6,ielem) = ip+(2*ny-1)

      if ( s_eqi ( region,'cavity')) then

        isotri(ielem) = 0

      else if ( s_eqi ( region,'cavity2')) then

        isotri(ielem) = 0

      else if ( s_eqi ( region,'channel')) then

        if ( ibump == 0) then

          if ( nbleft <= icol .and. icol < nbrite ) then
            isotri(ielem) = 1
          else
            isotri(ielem) = 0
          end if

        else if ( ibump == 1) then

          if ( irow == 1 .and. nbleft <= icol .and. icol < nbrite ) then
            isotri(ielem) = 2
          else if ( nbleft <= icol .and. icol < nbrite ) then
            isotri(ielem) = 1
          else
            isotri(ielem) = 0
          end if

        else if ( ibump == 2) then

          if ( nbleft <= icol .and. icol < nbrite ) then
            isotri(ielem) = 2
          else
            isotri(ielem) = 0
          end if

        else

          isotri(ielem) = 2

        end if

      else if ( s_eqi ( region,'step')) then

        isotri(ielem) = 0

      end if

    end if

    if ( maxnfl < neqnfl+2 ) then
      write ( *, * ) ' '
      write ( *, * ) 'SetNod - Fatal error!'
      write ( *, * ) '  Too many unknowns!'
      write ( *, * ) '  Processing node IP = ',ip
      write ( *, * ) '  The maximum allowed is MAXNFL = ',maxnfl
      write ( *, * ) '  This problem requires NEQNFL = ',neqnfl+2
      stop
    end if
!
!  Now determine what equations to associate with this node.
!
!  CAVITY:
!
    if ( s_eqi ( region,'cavity')) then
!
!  Top row:
!
      if ( irow == 2*ny-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UI'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VI'
!
!  Left, right, or bottom:
!
      else if ( icol == 1.or.icol.eq.2*nx-1.or.irow.eq.1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UW'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VW'
!
!  Interior:
!
      else

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'U'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'V'

      end if
!
!  CAVITY2:
!
    else if ( s_eqi ( region,'cavity2')) then
!
!  Top or bottome row:
!
      if ( irow == 1.or.irow.eq.2*ny-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UI'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VI'
!
!  Left or right side:
!
      else if ( icol == 1.or.icol.eq.2*nx-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UW'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VW'
!
!  Interior:
!
      else

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'U'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'V'

      end if
!
!  CHANNEL:
!
    else if ( s_eqi ( region,'channel')) then
!
!  The node lies on the left hand inflow boundary.
!  The horizontal and vertical velocities are specified.
!
      if ( icol == 1 .and. 1 < irow .and. irow < 2*ny-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UI'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VI'
!
!  The node lies on the right hand boundary.
!  The horizontal velocity is an unknown, the vertical velocity is zero.
!
      else if ( icol == 2*nx-1.and.1 < irow.and.irow < 2*ny-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'U'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VW'
!
!  The node lies on the moving bump surface.
!  The horizontal and vertical velocities are zero.
!
      else if ( irow == 1 .and. nbleft < icol .and. icol < nbrite ) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UB'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VB'
!
!  The node lies on a fixed wall.
!  The horizontal and vertical velocities are zero.
!
      else if ( icol == 1 .or. icol .eq. 2*nx-1 .or. &
           (irow == 1 .and. icol <= nbleft).or. &
           (irow == 1 .and. nbrite <= icol ) .or. irow .eq. 2*ny-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UW'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VW'
!
!  The node is a normal interior node.
!  The horizontal and vertical velocities are unknown.
!
      else

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'U'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'V'

      end if
!
!  STEP:
!
!  The node lies on the left hand inflow boundary.
!  The horizontal and vertical velocities are specified.
!
    else if ( s_eqi ( region,'step')) then

      if ( icol == 1.and.1 < irow.and.irow < 2*ny-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UI'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VI'
!
!  The node lies on the right hand boundary, above the JROW row.
!  The horizontal velocity is an unknown, the vertical velocity is zero.
!
      else if ( icol == 2*nx-1.and.jrow < irow.and.irow < 2*ny-1) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'U'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VW'
!
!  The node lies on a fixed wall or step.
!  The horizontal and vertical velocities are zero.
!
      else if ( &
           (irow == 1 .and. icol <= jcol ) .or. &
           (irow <= jrow .and. icol == jcol ) .or. &
           (irow == jrow .and. jcol <= icol ).or. &
           (irow == 2*ny-1) ) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'UW'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'VW'
!
!  The node lies in the "dead" zone.
!
      else if ( irow < jrow .and. jcol < icol ) then

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'U0'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'V0'
!
!  The node is a normal interior node.
!  The horizontal and vertical velocities are unknown.
!
      else

        neqnfl = neqnfl+1
        indx(1,ip) = neqnfl
        eqn(neqnfl) = 'U'

        neqnfl = neqnfl+1
        indx(2,ip) = neqnfl
        eqn(neqnfl) = 'V'

      end if
    end if
!
!  On nodes in an odd row and column, add a pressure equation.
!
    if ( irow2 == 1.and.icol2 == 1) then

      neqnfl = neqnfl+1

      if ( maxnfl < neqnfl ) then
        write ( *, * ) ' '
        write ( *, * ) 'SetNod - Fatal error!'
        write ( *, * ) '  Too many unknowns!'
        write ( *, * ) '  Processing node IP = ',ip
        write ( *, * ) '  The maximum allowed is MAXNFL = ',maxnfl
        write ( *, * ) '  This problem requires NEQNFL = ',neqnfl
        stop
      end if

      indx(3,ip) = neqnfl

      eqn(neqnfl) = 'P'

      if ( s_eqi ( region,'step')) then
        if ( irow < jrow .and. jcol < icol ) then
          eqn(neqnfl) = 'P0'
        end if
      end if

    else
      indx(3,ip) = 0
    end if

  end do
!
!  The last equation, which is guaranteed to be a pressure equation,
!  is replaced by a pressure boundary condition, associated with
!  an unknown.  (Even though we know this pressure will be zero).
!
  eqn(neqnfl) = 'PB'

  return
end
subroutine setpfl(area,etaq,isotri,nelem,node,np,phifl,xc,xquad,xsiq,yc,yquad)
!
!*****************************************************************************80
!
!! SETPFL evaluates and saves the basis functions at all quadrature points.  
!
!
!  Discussion:
!
!    The basis functions are precomputed and saved in this way for efficiency.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  AREA   Input, real ( kind = 8 ) AREA(3,NELEM).
!
!         AREA contains the area of each element.  These values are
!         needed when computed the integrals associated with the
!         finite element method.
!
!         For runs in which the region is allowed to change from
!         step to step, AREA must be recalculated at each step.
!
!  ETAQ   Input, real ( kind = 8 ) ETAQ(3).
!         The "Eta" coordinates of the quadrature points.
!
!  ISOTRI Input, integer ( kind = 4 ) ISOTRI(NELEM).
!
!         0, the element is NOT isoparametric.  The six node
!         triangle has straight sides.
!
!         1, the element is isoparametric.  The six node triangle
!         has curved sides.  Many computations involving such an
!         element must be computed by using a reference triangle,
!         and evaluating the jacobian of a transformation between
!         that triangle and the element.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM), contains the numbers
!         of the nodes that make up each element.  Element number
!         I is associated with nodes NODE(1,I) through NODE(6,I).
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  PHIFL  Output, real ( kind = 8 ) PHIFL(3,6,10,NELEM). 
!
!         PHIFL contains the value of a finite element basis function, its
!         derivative, or other information, evaluated at the quadrature
!         points.
!
!         The meaning of the entry PHIFL(I,J,K,L) is as follows. 
!         For the quadrature point I, and basis function J, in element L,
!         PHIFL(I,J,K,L) represents the value of:
!
!           K =  1, W, the finite element basis function for velocities;
!           K =  2, dWdX, the X derivative of W;
!           K =  3, dWdY, the Y derivative of W;
!           K =  4, Q, the finite element basis function for pressures;
!           K =  5, dQdX, the X derivative of Q;
!           K =  6, dQdY, the Y derivative of Q;
!           K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!           K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!           K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!           K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!         In particular, PHIFL(I,J,K,L) is the value of the quadratic
!         basis function W associated with local node J in element L,
!         evaluated at quadrature point I.
!
!         Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!         since there are only three linear basis functions.
!
!  XC     Input, real ( kind = 8 ) XC(NP), contains the X coordinates
!         of the nodes.
!
!  XQUAD  Input, real ( kind = 8 ) XQUAD(3,NELEM), contains the
!         X coordinates  of the quadrature points in a given element.
!
!  XSIQ   Input, real ( kind = 8 ) XSIQ(3).
!         The "Xsi" coordinates of the quadrature points.
!
!  YC     Input, real ( kind = 8 ) YC(NP), contains the Y coordinates
!         of the nodes.
!
!  YQUAD  Input, real ( kind = 8 ) YQUAD(3,NELEM), contains the
!         Y coordinates of the quadrature points in a given element.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) det
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  real ( kind = 8 ) etaq(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) q
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(3,nelem)
  real ( kind = 8 ) xq
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(3,nelem)
  real ( kind = 8 ) yq
!
!  Consider a particular element,
!  and a particular quadrature point (XQ,YQ) in that element.
!
!  Compute, at (XQ,YQ), the local values of the jacobian matrix
!  and its determinant.
!
!  Adjust the AREA array  
!
  do ielem = 1, nelem

    do i = 1, 3

      xq = xquad(i,ielem)
      yq = yquad(i,ielem)

      if ( isotri(ielem) == 2) then
        eta = etaq(i)
        xsi = xsiq(i)
        call trans(det,detadx,detady,dxsidx,dxsidy,eta,ielem, &
          nelem,node,np,xc,xsi,yc)
        area(i,ielem) = det*area(i,ielem)
      else
        eta = 0.0D+00
        xsi = 0.0D+00
      end if
!
!  Now consider each of the basis functions associated with a
!  node in the given element.
!
      do iq = 1,6
!
!  If the element is NOT isoparametric, compute the basis values
!  directly.
!
!  For isoparametric elements, use the reference triangle method.
!
        if ( isotri(ielem) == 0.or.isotri(ielem) == 1) then

          call bsp(q,dqdx,dqdy,ielem,iq,nelem,node,np,xc,xq,yc,yq)

          call qbf(ielem,iq,w,dwdx,dwdy,nelem,node,np,xc,xq,yc,yq)

          dxsidx = 1.0D+00
          dxsidy = 0.0D+00
          detadx = 0.0D+00
          detady = 1.0D+00

        else

          call refqbf(w,dwdx,dwdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)

          call refbsp(q,dqdx,dqdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)

        end if
!
!  Store the values into PHIFL.
!
        phifl(i,iq,1,ielem) = w
        phifl(i,iq,2,ielem) = dwdx
        phifl(i,iq,3,ielem) = dwdy
        phifl(i,iq,4,ielem) = q
        phifl(i,iq,5,ielem) = dqdx
        phifl(i,iq,6,ielem) = dqdy
        phifl(i,iq,7,ielem) = dxsidx
        phifl(i,iq,8,ielem) = dxsidy
        phifl(i,iq,9,ielem) = detadx
        phifl(i,iq,10,ielem) = detady    

      end do
    end do
  end do

  return
end
subroutine setprb(eqn,indx,maxcofrb,maxelm,maxnfl,nelem,neqnfl, &
  ncofrb,node,np,phifl,phirb,rb)
!
!*****************************************************************************80
!
!! SETPRB evaluates the reduced basis functions at each quadrature point.
!
!
!  Discussion:
!
!    The routine is given:
!
!      GFLRB, the full solution at which the reduced basis was generated;
!      PHIFL, the value of the finite element basis functions
!        at each quadrature point,
!      RB, the reduced basis vectors;
!
!    and computes:
!
!      PHIRB, the value of the reduced basis functions at each quadrature
!        point, for each reduced basis vector.
!
!    Note that the PHIFL contains the values of the finite element basis
!    functions at each quadrature point, and so we can compute ANY possible
!    finite element solution by "multiplying" PHIFL by a choice of coefficients
!    GFL.
!
!    What we are essentially doing is picking particular choices of coefficients,
!    namely the columns of the RB array, and computing the resulting values
!    of U, V, and P at each quadrature point in each element.  Then, later
!    on, a linear combination of reduced basis vectors can be evaluated
!    easily at any quadrature point simply by "multiplying" the reduced
!    basis coefficients by the entries of PHIRB that tell us the values
!    of the basis functions associated with each reduced basis vector.
!
!    Note that although the finite element basis functions are the same
!    for the U and V velocities, this is NOT true in the reduced basis.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the global index of U,
!         V and P at that node, or 0 or a negative value.  The global
!         index of U, V, or P is the index of the coefficient vector
!         that contains the value of the finite element coefficient
!         associated with the corresponding basis function at the
!         given node.
! 
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  MAXNFL Input, integer ( kind = 4 ) MAXNFL.
!         MAXNFL is the maximum number of equations or coefficients allowed
!         for the full system.  MAXNFL must be used instead of NEQNFL as
!         the leading dimension of certain multi-dimensional arrays.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB.
!         NCOFRB is the number of equations in the reduced system,
!         and also the number of reduced basis coefficients
!         which need to be determined by those equations.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!
!         NODE(I,J) contains, for an element J, the global index of
!         the node whose local number in J is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!           Global nodes   Elements      NODE
!                                                        1  2  3  4  5  6
!           74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                          |  /   /|
!           73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                          |/   /  |
!           72  82  92     2   1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP.
!         NP is the number of nodes used to define the finite element mesh.
!         Typically, the mesh is generated as a rectangular array, with
!         an odd number of nodes in the horizontal and vertical directions.
!         The formula for NP is NP = (2*NX-1)*(2*NY-1).
! 
!  PHIFL  Input, real ( kind = 8 ) PHIFL(3,6,10,NELEM). 
!
!         PHIFL contains the value of a finite element basis function, its
!         derivative, or other information, evaluated at the quadrature
!         points.
!
!         The meaning of the entry PHIFL(I,J,K,L) is as follows. 
!         For the quadrature point I, and basis function J, in element L,
!         PHIFL(I,J,K,L) represents the value of:
!
!           K =  1, W, the finite element basis function for velocities;
!           K =  2, dWdX, the X derivative of W;
!           K =  3, dWdY, the Y derivative of W;
!           K =  4, Q, the finite element basis function for pressures;
!           K =  5, dQdX, the X derivative of Q;
!           K =  6, dQdY, the Y derivative of Q;
!           K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!           K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!           K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!           K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!         In particular, PHIFL(I,J,K,L) is the value of the quadratic
!         basis function W associated with local node J in element L,
!         evaluated at quadrature point I.
!
!         Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!         since there are only three linear basis functions.
!
!  PHIRB  Output, real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!  RB     Input, real ( kind = 8 ) rb(maxnfl,ncofrb).
!         RB is the NEQNFL by NCOFRB array of reduced basis vectors.
!
!         RB is generated by computing a finite element solution GFL,
!         which is saved for later reference as "GFLRB".
!         GFLRB is copied into the first column of RB.
!         Then, we compute the first NCOFRB derivatives of GFLRB with
!         respect to a parameter.  The first derivative
!         is stored in column 1 of RB, and so on.  Then we orthogonalize
!         the columns of RB.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dqrbdx
  real ( kind = 8 ) dqrbdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) dwu0rbdx
  real ( kind = 8 ) dwu0rbdy
  real ( kind = 8 ) dwurbdx
  real ( kind = 8 ) dwurbdy
  real ( kind = 8 ) dwv0rbdx
  real ( kind = 8 ) dwv0rbdy
  real ( kind = 8 ) dwvrbdx
  real ( kind = 8 ) dwvrbdy
  character ( len = 2 ) eqn(neqnfl)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ieqnrb
  integer ( kind = 4 ) iglob
  integer ( kind = 4 ) ilocal
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) iquad
  logical s_eqi
  integer ( kind = 4 ) nglob
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) q
  real ( kind = 8 ) qrb
  real ( kind = 8 ) rb(maxnfl,ncofrb)
  real ( kind = 8 ) w
  real ( kind = 8 ) wu0rb
  real ( kind = 8 ) wurb
  real ( kind = 8 ) wv0rb
  real ( kind = 8 ) wvrb
!
!  Consider element IELEM...
!
  do ielem = 1,nelem
!
!  ...and quadrature point IQUAD in element IELEM...
!
    do iquad = 1,3
!
!  ...and reduced basis function IEQNRB, evaluated at
!  quadrature point IQUAD in element IELEM.
!
      do ieqnrb = 1,ncofrb

        wurb = 0.0D+00
        dwurbdx = 0.0D+00
        dwurbdy = 0.0D+00

        wvrb = 0.0D+00
        dwvrbdx = 0.0D+00
        dwvrbdy = 0.0D+00

        qrb = 0.0D+00
        dqrbdx = 0.0D+00
        dqrbdy = 0.0D+00

        wu0rb = 0.0D+00
        dwu0rbdx = 0.0D+00
        dwu0rbdy = 0.0D+00

        wv0rb = 0.0D+00
        dwv0rbdx = 0.0D+00
        dwv0rbdy = 0.0D+00
!
!  Now add up the U, V, or P finite element coefficients, weighted by the
!  values of the finite element basis functions or derivatives, at the
!  quadrature point.
!
        do ilocal = 1,6

          w    = phifl(iquad,ilocal,1,ielem)
          dwdx = phifl(iquad,ilocal,2,ielem)
          dwdy = phifl(iquad,ilocal,3,ielem)

          q    = phifl(iquad,ilocal,4,ielem)
          dqdx = phifl(iquad,ilocal,5,ielem)
          dqdy = phifl(iquad,ilocal,6,ielem)

          nglob = node(ilocal,ielem)
         
          iglob = indx(1,nglob)

          wurb    = wurb   +rb(iglob,ieqnrb)*w
          dwurbdx = dwurbdx+rb(iglob,ieqnrb)*dwdx
          dwurbdy = dwurbdy+rb(iglob,ieqnrb)*dwdy

          if ( s_eqi ( eqn(iglob),'u')) then
            wu0rb    = wu0rb   +rb(iglob,ieqnrb)*w
            dwu0rbdx = dwu0rbdx+rb(iglob,ieqnrb)*dwdx
            dwu0rbdy = dwu0rbdy+rb(iglob,ieqnrb)*dwdy
          end if

          iglob = indx(2,nglob)

          wvrb    = wvrb   +rb(iglob,ieqnrb)*w
          dwvrbdx = dwvrbdx+rb(iglob,ieqnrb)*dwdx
          dwvrbdy = dwvrbdy+rb(iglob,ieqnrb)*dwdy

          if ( s_eqi ( eqn(iglob),'v')) then
            wv0rb    = wv0rb   +rb(iglob,ieqnrb)*w
            dwv0rbdx = dwv0rbdx+rb(iglob,ieqnrb)*dwdx
            dwv0rbdy = dwv0rbdy+rb(iglob,ieqnrb)*dwdy
          end if

          iglob = indx(3,nglob)

          if ( 0 < iglob ) then

            qrb    = qrb   +rb(iglob,ieqnrb)*q
            dqrbdx = dqrbdx+rb(iglob,ieqnrb)*dqdx
            dqrbdy = dqrbdy+rb(iglob,ieqnrb)*dqdy

          end if

        end do
!
!  Save the values of the finite element basis functions associated
!  with the given reduced basis vector.
!
        phirb(iquad,ieqnrb,1,ielem) = wurb
        phirb(iquad,ieqnrb,2,ielem) = dwurbdx
        phirb(iquad,ieqnrb,3,ielem) = dwurbdy 

        phirb(iquad,ieqnrb,4,ielem) = wvrb
        phirb(iquad,ieqnrb,5,ielem) = dwvrbdx
        phirb(iquad,ieqnrb,6,ielem) = dwvrbdy

        phirb(iquad,ieqnrb,7,ielem) = qrb
        phirb(iquad,ieqnrb,8,ielem) = dqrbdx
        phirb(iquad,ieqnrb,9,ielem) = dqrbdy

        phirb(iquad,ieqnrb,10,ielem) = wu0rb
        phirb(iquad,ieqnrb,11,ielem) = dwu0rbdx
        phirb(iquad,ieqnrb,12,ielem) = dwu0rbdy 

        phirb(iquad,ieqnrb,13,ielem) = wv0rb
        phirb(iquad,ieqnrb,14,ielem) = dwv0rbdx
        phirb(iquad,ieqnrb,15,ielem) = dwv0rbdy
      end do
    end do

  end do

  return
end
subroutine setq3(area,etaq,isotri,nelem,node,np,wquad,xc,xquad,xsiq,yc,yquad)
!
!*****************************************************************************80
!
!! SETQ3 sets the abscissas and weights for a three point quadrature
!  rule on a triangle.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  AREA   Output, real ( kind = 8 ) AREA(3,NELEM).
!
!  ETAQ   Input, real ( kind = 8 ) ETAQ(3).
!         The "Eta" coordinates of the quadrature points.
!
!  ISOTRI Input, integer ( kind = 4 ) ISOTRI(NELEM).
!
!         0, the element is NOT isoparametric, and the nodes never move.
!         That means that the quadrature points are only computed once.
!
!         1, the element is NOT isoparametric, but the nodes may move.
!         Quadrature point locations must be updated on each step.
!         This could occur for elements above, but not touching, the bump.
!
!         2, the element is isoparametric.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM).
!
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  WQUAD  Input, real ( kind = 8 ) WQUAD(3), the weights for Gaussian
!         quadrature.
!
!  XC     Input, real ( kind = 8 ) XC(NP).
!
!         The X coordinates of the nodes.
!
!  XQUAD  Output, real ( kind = 8 ) XQUAD(3,NELEM).
!
!         The X coordinates of the quadrature points for
!         each element.
!
!  XSIQ   Input, real ( kind = 8 ) XSIQ(3).
!         The "Xsi" coordinates of the quadrature points.
!
!  YC     Input, real ( kind = 8 ) YC(NP).
!
!         The Y coordinates of the nodes.
!
!  YQUAD  Output, real ( kind = 8 ) YQUAD(3,NELEM).
!
!         The Y coordinates of the quadrature points for
!         each element.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) area(3,nelem)
  real ( kind = 8 ) eta
  real ( kind = 8 ) etaq(3)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ip1
  integer ( kind = 4 ) ip2
  integer ( kind = 4 ) ip3
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) wquad(3)
  real ( kind = 8 ) x
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xquad(3,nelem)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) xsiq(3)
  real ( kind = 8 ) y
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yquad(3,nelem)
!
!  Set the weights.
!
  do i = 1,3
    wquad(i) = 1.0D+00 / 6.0D+00
  end do
!
!  Set the quadrature points in the reference element.
!
  xsiq(1) = 0.5D+00
  etaq(1) = 0.5D+00

  xsiq(2) = 1.0D+00
  etaq(2) = 0.5D+00

  xsiq(3) = 0.5D+00
  etaq(3) = 0.0D+00
!
!  Set the X, Y coordinates of quadrature points for each element.
!
  do ielem = 1,nelem

    do i = 1,3
      xsi = xsiq(i)
      eta = etaq(i)
      call xofxsi(eta,ielem,nelem,node,np,x,xc,xsi,y,yc)
      xquad(i,ielem) = x
      yquad(i,ielem) = y
    end do
!
!  We only calculate true areas for nonisoparametric elements.
!
    ip1 = node(1,ielem)
    ip2 = node(2,ielem)
    ip3 = node(3,ielem)

    do iquad = 1,3

      if ( isotri(ielem) == 0.or.isotri(ielem).eq.1) then

        area(iquad,ielem) = wquad(iquad)*abs( &
              (yc(ip1)+yc(ip2))*(xc(ip2)-xc(ip1)) &
             +(yc(ip2)+yc(ip3))*(xc(ip3)-xc(ip2)) &
             +(yc(ip3)+yc(ip1))*(xc(ip1)-xc(ip3)) )

      else

        area(iquad,ielem) = wquad(iquad)

      end if
 
    end do

  end do

  return
end
subroutine setxy(gridx,gridy,ibs,np,nparb,nx,ny,region,splbmp, &
  taubmp,xbl,xbr,xc,xrange,ybl,ybr,yc,yrange)
!
!*****************************************************************************80
!
!! SETXY sets the X and Y coordinates of the nodes.
!
!
!  SETXY assumes that the nodes are numbered
!  in "stacks", starting with the least X and Y coordinates,
!  then fixing X and running through all values of Y, then
!  increasing X to the next value and running through all
!  values of Y, and so on.  For example:
!
!    5  10  15
!    4   9  14
!    3   8  13
!    2   7  12
!    1   6  11
!
!  SETXY allows a certain number of schemes for computing the
!  grid in the X and Y directions, aside from uniform spacing.
!  However, SETXY forces the nodes in "even" rows to have the
!  Y coordinates that are the average of the nodes above and
!  below, and nodes in "even" columns to have the X coordinates
!  that are the average of the nodes left and right of them.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GRIDX 
!    Input, character ( len = 20 ) GRIDX.
!    GRIDX tells how the finite element nodes should be layed out
!    in the X direction.
!    'uniform' makes them equally spaced.
!    'cos' uses the COS function to cluster them near edges.
!    'sqrtsin' uses the SQRT(SIN()) function.
!
!  GRIDY 
!    Input, character ( len = 20 ) GRIDY.
!    GRIDY tells how the finite element nodes should be layed out
!    in the Y direction.
!    'uniform' makes them equally spaced.
!    'cos' uses the COS function to cluster them near edges.
!    'sqrtsin' uses the SQRT(SIN()) function.
!
!  IBS   
!    Input, integer ( kind = 4 ) IBS.
!    IBS is the bump shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!  NP   
!    Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  NPARB 
!    Input, integer ( kind = 4 ) NPARB.
!
!    The number of parameters associated with the position and
!    shape of the bump.
!
!    Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NX   
!    Input, integer ( kind = 4 ) NX.
!
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!
!    Roughly speaking, NX (or 2*NX) is the number of elements along
!    a line in the X direction.
!
!  NY    
!    Input, integer ( kind = 4 ) NY.
!
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!
!    Roughly speaking, NY (or 2*NY) is the number of elements along
!    a line in the Y direction.
!
!  REGION
!    Input, character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'cavity2', a driven cavity, 1 unit on each side, open on
!    the top and bottom with tangential velocity specifications there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!  SPLBMP
!    Input, real ( kind = 8 ) SPLBMP(NPARB+2).
!
!    SPLBMP contains the spline coefficients for the bump.
!
!  TAUBMP
!    Input, real ( kind = 8 ) TAUBMP(NPARB+2).
!
!    TAUBMP contains the location of the spline abscissas for
!    the bump.  There are NPARB+2 of them, because the end values
!    of the spline are constrained to have particular values.
!
!  XBL  
!    Input, real ( kind = 8 ) XBL.
!    XBL is the X coordinate of the left corner
!    of the bump.
!
!  XBR  
!    Input, real ( kind = 8 ) XBR.
!    XBR is the X coordinate of the right corner
!    of the bump.
!
!  XC   
!    Output, real ( kind = 8 ) XC(NP).
!    XC is the X coordinates of the nodes.
!
!  YBL   
!    Input, real ( kind = 8 ) YBL.
!    YBL is the Y coordinate of the left corner of the bump.
!
!  YBR   
!    Input, real ( kind = 8 ) YBR.
!    YBR is the Y coordinate of the right corner of the bump.
!
!  YC    
!    Output, real ( kind = 8 ) YC(NP).
!    YC is the Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
!
  character ( len = 20 ) gridx
  character ( len = 20 ) gridy
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) icol
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) ipe
  integer ( kind = 4 ) ipn
  integer ( kind = 4 ) ips
  integer ( kind = 4 ) ipw
  integer ( kind = 4 ) irow
  integer ( kind = 4 ) jcol
  integer ( kind = 4 ) jrow
  logical s_eqi
  character ( len = 20 ) region
  real ( kind = 8 ) splbmp(nparb+2)
  real ( kind = 8 ) taubmp(nparb+2)
  real ( kind = 8 ) x
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xhi
  real ( kind = 8 ) xlo
  real ( kind = 8 ) xrange
  real ( kind = 8 ) y
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybot
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yhi
  real ( kind = 8 ) ylo
  real ( kind = 8 ) yrange
!
  ip = 0
!
!  For the step problem only,
!  set the odd row and column of the corner of the step.
!
  if ( s_eqi ( region,'step')) then
    jcol = 2*nint(xbr*(nx-1)/xrange)+1
    jrow = 2*nint(ybr*(ny-1)/yrange)+1
  end if
!
!  Consider each column of the region.
!
  do icol = 1,2*nx-1

    if ( s_eqi ( region,'cavity').or.s_eqi ( region,'cavity2')) then

      ilo = 1
      ihi = 2*nx-1
      xlo = 0.0D+00
      xhi = xrange

      call grid(gridx,icol,ihi,ilo,x,xhi,xlo)

    else if ( s_eqi ( region,'channel')) then

      ilo = 1
      ihi = 2*nx-1
      xlo = 0.0D+00
      xhi = xrange

      call grid(gridx,icol,ihi,ilo,x,xhi,xlo)

      if ( abs(x-xbl)*(2*nx-2) <= 0.5) then
        x = xbl
      else if ( abs(x-xbr)*(2*nx-2) <= 0.5) then
        x = xbr
      end if

    else if ( s_eqi ( region,'step')) then

      if ( icol < jcol) then

        ilo = 1
        ihi = jcol
        xlo = 0.0D+00
        xhi = xbr

        call grid(gridx,icol,ihi,ilo,x,xhi,xlo)

      else if ( icol == jcol) then
        x = xbr
      else if ( jcol < icol ) then

        ilo = jcol
        ihi = 2*nx-1
        xlo = xbr
        xhi = xrange

        call grid(gridx,icol,ihi,ilo,x,xhi,xlo)

      end if

    end if
!
!  Consider each row of the region.
!
    do irow = 1,2*ny-1

      ip = ip+1

      if ( s_eqi ( region,'cavity').or.s_eqi ( region,'cavity2')) then

        ilo = 1
        ihi = 2*ny-1
        ylo = 0.0D+00
        yhi = yrange

        call grid(gridy,irow,ihi,ilo,y,yhi,ylo)

      else if ( s_eqi ( region,'channel')) then

        if ( x <= xbl) then
          ybot = ybl
        else if ( xbl <= x .and. x <= xbr ) then
          if ( ibs == 0) then
            call pcval(nparb+1,x,taubmp,ybot,splbmp)
          else if ( ibs == 1) then
            call plval(nparb+2,x,taubmp,ybot,splbmp)
          else if ( ibs == 2) then
            call pqval(nparb+2,x,taubmp,ybot,splbmp)
          end if
        else
          ybot = ybr
        end if

        ilo = 1
        ihi = 2*ny-1
        ylo = ybot
        yhi = yrange

        call grid(gridy,irow,ihi,ilo,y,yhi,ylo)

      else if ( s_eqi ( region,'step')) then

        if ( irow < jrow) then

          ilo = 1
          ihi = jrow
          ylo = 0.0D+00
          yhi = ybr

          call grid(gridy,irow,ihi,ilo,y,yhi,ylo)

        else if ( irow == jrow) then
          y = ybr
        else if ( jrow < irow ) then

          ilo = jrow
          ihi = 2*ny-1
          ylo = ybr
          yhi = yrange

          call grid(gridy,irow,ihi,ilo,y,yhi,ylo)

        end if

      end if

      xc(ip) = x
      yc(ip) = y

    end do
  end do
!
!  Average the X coordinates of all nodes that lie in even columns.
!
  do irow = 1,2*ny-1
    do icol = 2,2*nx-2,2
      ip = (icol-1)*(2*ny-1)+irow
      ipw = ip-(2*ny-1)
      ipe = ip+(2*ny-1)
      xc(ip) = 0.5D+00 *(xc(ipe)+xc(ipw))
    end do
  end do
!
!  Average the Y coordinates of all nodes that lie in even rows.
!
  do irow = 2,2*ny-2,2
    do icol = 1,2*nx-1
      ip = (icol-1)*(2*ny-1)+irow
      ipn = ip+1
      ips = ip-1
      yc(ip) = 0.5D+00 *(yc(ipn)+yc(ips))
    end do
  end do

  return
end
subroutine step(ibs,ibump,ifs,iopt,maxopt,maxpar,nbcrb,npar,nparb,nparf, &
  npe,nx,ny,par,region,reynld,tolnew,tolopt,tolsim,wateb,watep,wateu, &
  watev,xbl,xbr,xprof,xrange,ybl,ybr,yrange)
!
!*****************************************************************************80
!
!! STEP sets up a forward facing step problem.
!
!
!  Reference:
!
!    Janet Peterson,
!    The Reduced Basis Method for Incompressible Viscous Flow Calculations,
!    SIAM Journal of Scientific and Statistical Computing,
!    Volume 10, Number 4, pages 777-786, July 1989.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  IBS   
!    integer ( kind = 4 ) IBS.
!    IBS is the bump shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!  IBUMP 
!    integer ( kind = 4 ) IBUMP.
!    IBUMP determines where isoparametric elements will be used.
!
!    0, no isoparametric elements will be used.
!       The Y coordinates of midside nodes of elements above the
!       bump will be recomputed so that the sides are straight.
!
!    1, isoparametric elements will be used only for the
!       elements which directly impinge on the bump.
!       Midside nodes of nonisoparametric elements above the
!       bump will be recomputed so that the sides are straight.
!
!    2, isoparametric elements will be used for all elements
!       which are above the bump.  All nodes above the bump
!       will be equally spaced in the Y direction.
!
!    3, isoparametric elements will be used for all elements.
!       All nodes above the bump will be equally spaced in
!       the Y direction.
!
!  IFS   
!    integer ( kind = 4 ) IFS.
!    IFS is the inflow shape option.
!    0, piecewise constant function.
!    1, piecewise linear function.
!    2, piecewise quadratic function.
!
!  IOPT  
!    integer ( kind = 4 ) IOPT(MAXPAR).
!    IOPT is used during an optimization.  For each parameter I,
!    the meaning of IOPT(I) is:
!    0, the parameter value must remain fixed;
!    1, the parameter value may be varied.
!
!  MAXOPT
!    integer ( kind = 4 ) MAXOPT.
!    MAXOPT is the maximum number of optimization steps.
!
!  MAXPAR
!    integer ( kind = 4 ) MAXPAR.
!    MAXPAR is the maximum number of parameters allowed.
!    MAXPAR = MAXPARF + MAXPARB + 1.
!
!  NBCRB 
!    integer ( kind = 4 ) NBCRB.
!    NBCRB is the number of independent boundary condition
!    vectors used for the reduced basis.  NBCRB is normally
!    at least 1, and must be no more than MAXBCRB.
!
!  NPAR  
!    integer ( kind = 4 ) NPAR.
!    NPAR is the number of parameters.
!      NPAR = NPARF + NPARB + 1.
!    The parameters control the shape and strength of the inflow,
!    the shape of the bump, and the value of the Reynolds number.
!
!  NPARB 
!    integer ( kind = 4 ) NPARB.
!    NPARB is the number of parameters associated with the position and
!    shape of the bump.
!
!    Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF 
!    integer ( kind = 4 ) NPARF.
!    NPARF is the number of parameters associated with the
!    inflow.  NPARF must be at least 1.
!
!  NPE
!    integer ( kind = 4 ) NPE.
!    NPE is the number of nodes per element.
!
!  NX    
!    integer ( kind = 4 ) NX.
!    NX controls the spacing of nodes and elements in
!    the X direction.  There are 2*NX-1 nodes along various
!    lines in the X direction.
!
!    The number of elements along a line in the X direction is
!    NX-1 (or 2*(NX-1) to make a full rectangular strip).
!
!  NY    
!    integer ( kind = 4 ) NY.
!    NY controls the spacing of nodes and elements in
!    the Y direction.  There are 2*NY-1 nodes along various
!    lines in the Y direction.
!
!    The number of elements along a line in the Y direction is
!    NY-1 (or 2*(NY-1) to make a full vertical strip).
!
!  PAR   
!    real ( kind = 8 ) PAR(NPAR).
!    PAR contains the values of the problem parameters.
!
!      PAR(1:NPARF)             = inflow controls.
!      PAR(NPARF+1:NPARF+NPARB) = bump controls.
!      PAR(NPARF+NPARB+1)       = the REYNLD parameter.
!
!  REGION
!    character ( len = 20 ) REGION.
!    REGION specifies the flow region.
!
!    'cavity', a driven cavity, 1 unit on each side, open on
!    the top with a tangential velocity specification there.
!
!    'channel', a channel, 10 units long by 3 high, inflow on
!    the left, outflow on the right, with a bump on the bottom.
!
!    'step', a channel, 12 units long by 3 high, inflow on the
!    left, outflow on the right, with a step on the bottom.
!
!  REYNLD
!    real ( kind = 8 ) REYNLD.
!    REYNLD is the current value of the Reynolds number.
!    Normally, REYNLD is stored as PARA(NPARF+NPARB+1).
!
!  TOLNEW
!    real ( kind = 8 ) TOLNEW.
!    TOLNEW is the convergence tolerance for the Newton iteration.
!
!  TOLOPT
!    real ( kind = 8 ) TOLOPT.
!    TOLOPT is the convergence tolerance for the optimization.
!
!  TOLSIM
!    real ( kind = 8 ) TOLSIM.
!    TOLSIM is the convergence tolerance for the Picard iteration.
!
!  WATEB 
!    real ( kind = 8 ) WATEB.
!    WATEB is the multiplier of the bump control cost used
!    when computing the total cost.
!
!  WATEP,
!  WATEU,
!  WATEV 
!    real ( kind = 8 ) WATEP, WATEU, WATEV.
!
!    WATEP, WATEU and WATEV are weights used in computing the
!    cost function based on the costs of the flow discrepancy.
!
!  XBL   
!    real ( kind = 8 ) XBL.
!    XBL is the X coordinate of the left corner of the bump.
!
!  XBR   
!    real ( kind = 8 ) XBR.
!    XBR is the X coordinate of the right corner of the bump.
!
!  XPROF 
!    real ( kind = 8 ) XPROF.
!    XPROF is the X coordinate at which the profile is measured. 
!    XPROF should be a grid value!
!
!  XRANGE
!    real ( kind = 8 ) XRANGE.
!    XRANGE is the total width of the region.
!
!  YBL   
!    real ( kind = 8 ) YBL.
!    YBL is the Y coordinate of the left corner of the bump.
!
!  YBR   
!    real ( kind = 8 ) YBR.
!    YBR is the Y coordinate of the right corner of the bump.
!
!  YRANGE
!    real ( kind = 8 ) YRANGE.
!    YRANGE is the total height of the region.
!
  implicit none
!
  integer ( kind = 4 ) maxpar
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibs
  integer ( kind = 4 ) ibump
  integer ( kind = 4 ) ifs
  integer ( kind = 4 ) iopt(maxpar)
  integer ( kind = 4 ) maxopt
  integer ( kind = 4 ) nbcrb
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nparf
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nx
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(maxpar)
  character ( len = 20 ) region
  real ( kind = 8 ) reynld
  real ( kind = 8 ) tolnew
  real ( kind = 8 ) tolopt
  real ( kind = 8 ) tolsim
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) xprof
  real ( kind = 8 ) xrange
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yrange
!
  ibs = 0
  ibump = 0
!
!  The inflow is to be a (piecewise) quadratic function.
!
  ifs = 2
  maxopt = 10
  nbcrb = 1
  nparb = 0
!
!  There is only one unknown coefficient for the inflow function,
!  so there is only one "piece".
!
  nparf = 1
  npe = 6
!
!  Peterson used a nonuniform mesh with NX = 70 and NY=35!
!
  nx = 11
  ny = 4
  region = 'step'
  tolnew = 0.0000000001D+00
  tolopt = 0.000000001D+00
  tolsim = 0.0000000001D+00
  wateb = 0.0D+00
  wateu = 1.0D+00
  watev = 1.0D+00
  watep = 0.0D+00
  xbl = 4.0D+00
  xbr = 4.0D+00
  xprof = 3.0D+00
  xrange = 12.0D+00
  ybl = 1.0D+00
  ybr = 1.0D+00
  yrange = 3.0D+00
!
!  Set things that depend on other things.
!
  npar = nparf+nparb+1

  do i = 1,nparf
    iopt(i) = 0
  end do

  do i = nparf+1,nparf+nparb
    iopt(i) = 0
  end do

  iopt(nparf+nparb+1) = 1
!
!  The inflow parameter should be 1 to match Peterson's paper.
!
  par(1) = 1.0D+00
!
!  The REYNLD parameter may be varied.  Here, it is arbitrarily
!  set to 1.  Peterson worked with values up to 1500.
!
  reynld = 1.0D+00
  par(2) = reynld

  return
end
subroutine target(cost0,gfl,gfltar,indx,maxnfl,maxny,maxparb, &
  neqnfl,np,npar,nparb,nprof,ny,par,partar,splbmp, &
  taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)
!
!*****************************************************************************80
!
!! TARGET is called to save the current parameters and solution
!  as the "target solution".
!
!
!  TARGET was modified on 08 August 1996.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 August 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  COST0  Output, real ( kind = 8 ) COST0.
!         COST0 is the cost of the solution with PAR = GFL=0.
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!         GFL contains the current solution estimate for the full problem,
!         containing the pressure and velocity coefficients.
!         The vector INDX must be used to index this data.
!
!  GFLTAR Output, real ( kind = 8 ) GFLTAR(NEQNFL).
!         GFLTAR is a target solution, used to generate data that defines
!         the cost functional.  The corresponding parameters are PARTAR.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!         INDX(I,J) contains, for each node J, the global index of U,
!         V and P at that node, or 0 or a negative value.  The global
!         index of U, V, or P is the index of the coefficient vector
!         that contains the value of the finite element coefficient
!         associated with the corresponding basis function at the
!         given node.
!
!  MAXNFL Input, integer ( kind = 4 ) MAXNFL, the maximum number of equations in the
!         full system.
!
!  MAXNY  Input, integer ( kind = 4 ) MAXNY.
!         MAXNY is the maximum size of NY that the program can handle.
!
!  MAXPARB
!         Input, integer ( kind = 4 ) MAXPARB.
!         The maximum number of bump parameters allowed.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of finite element equations used
!         to define the horizontal and vertical velocities and the
!         pressure.
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  NPAR   Input, integer ( kind = 4 ) NPAR.
!
!         The number of parameters.  NPAR = NPARF + NPARB + 1.
!
!         The parameters control the shape of the inflow,
!         the shape of the bump obstacle, and the strength of the
!         flow.
!
!  NPARB  Input, integer ( kind = 4 ) NPARB.
!
!         The number of parameters associated with the position and
!         shape of the bump.
!
!         Note that if NPARB = 0, the bump is replaced by a flat wall.
!
!  NPARF  Input, integer ( kind = 4 ) NPARF.
!
!         NPARF is the number of parameters associated with the
!         inflow.  NPARF must be at least 1. 
!
!  NY     Input, integer ( kind = 4 ) NY, the number of nodes in each column.
!
!  PAR    Input, real ( kind = 8 ) PAR(NPAR).
!
!         PAR is the current set of parameter values, including the
!         Reynolds parameter, the flow parameters, and the bump parameters.
!
!  PARTAR Output, real ( kind = 8 ) PARTAR(NPAR).
!         PARTAR is the value of the parameters that generated the
!         target solution contained in GFLTAR.
!
!  SPLBMP Input, real ( kind = 8 ) SPLBMP(NPARB+2).
!         SPLBMP contains the spline coefficients for the bump.
!
!  TAUBMP Input, real ( kind = 8 ) TAUBMP(NPARB+2).
!         TAUBMP contains the location of the spline abscissas for
!         the bump.  There are NPARB+2 of them, because the end values
!         of the spline are constrained to have particular values.
!
!  WATEB  Input, real ( kind = 8 ) WATEB.
!         WATEB is the multiplier of the bump control cost used
!         when computing the total cost.
!
!  WATEP,
!  WATEU,
!  WATEV  Input, real ( kind = 8 ) WATEP, WATEU, WATEV.
!         WATEP, WATEU and WATEV are weights used in computing the 
!         cost function based on the costs of the flow discrepancy.
!
!  XBL    Input, real ( kind = 8 ) XBL.
!         XBL is the X coordinate of the left corner of the bump.
!
!  XBR    Input, real ( kind = 8 ) XBR.
!         XBR is the X coordinate of the right corner of the bump.
!
!  YBL    Input, real ( kind = 8 ) YBL.
!         The Y coordinate of the left corner of the bump.
!
!  YBR    Input, real ( kind = 8 ) YBR.
!         YBR is the Y coordinate of the right corner of the bump.
!
!  YC     Input, real ( kind = 8 ) YC(NP).
!         YC contains the Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) maxny
  integer ( kind = 4 ) maxparb
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
!
!  Set parameters that are dependent on parameters that are dependent
!  on parameters.
!
  real ( kind = 8 ) cost
  real ( kind = 8 ) cost0
  real ( kind = 8 ) costb
  real ( kind = 8 ) costp
  real ( kind = 8 ) costu
  real ( kind = 8 ) costv
  real ( kind = 8 ) gfl(maxnfl)
  real ( kind = 8 ) gfltar(maxnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) nparb
  integer ( kind = 4 ) nprof(2*maxny-1)
  integer ( kind = 4 ) ny
  real ( kind = 8 ) par(npar)
  real ( kind = 8 ) partar(npar)
  real ( kind = 8 ) splbmp(maxparb+2)
  real ( kind = 8 ) taubmp(maxparb+2)
  real ( kind = 8 ) wateb
  real ( kind = 8 ) watep
  real ( kind = 8 ) wateu
  real ( kind = 8 ) watev
  real ( kind = 8 ) xbl
  real ( kind = 8 ) xbr
  real ( kind = 8 ) ybl
  real ( kind = 8 ) ybr
  real ( kind = 8 ) yc(np)
!
  do i = 1,npar
    partar(i) = 0.0D+00
  end do

  do i = 1,neqnfl
    gfltar(i) = 0.0D+00
  end do

  call getcst(cost,costb,costp,costu,costv,gfl,gfltar, &
    indx,neqnfl,np,nparb,nprof,ny,splbmp, &
    taubmp,wateb,watep,wateu,watev,xbl,xbr,ybl,ybr,yc)

  cost0 = cost
  write ( *, * ) ' '
  write ( *, * ) '"Cost" of target GFLTAR versus zero GFL:',cost0

  do i = 1,npar
    partar(i) = par(i)
  end do

  do i = 1,neqnfl
    gfltar(i) = gfl(i)
  end do

  return
end
subroutine trans(det,detadx,detady,dxsidx,dxsidy,eta,ielem, &
  nelem,node,np,xc,xsi,yc)
!
!*****************************************************************************80
!
!! TRANS calculates the biquadratic transformation which maps the
!  reference element in (XSI,ETA) space into a particular
!  isoparametric element in (X,Y) space.
!
!  We know everything about the isoparametric element once we
!  specify the location of its six nodes.
!
!  TRANS computes the entries of the jacobian of the transformation
!  and the determinant of the jacobian.  Essentially, the jacobian
!  records the relationship between derivatives with respect to XSI
!  and ETA and a point in the reference element, and derivatives
!  with respect to X and Y of the same function as defined in the
!  isoparametric element.
!
!  The four entries of the jacobian are symbolically named DETADX,
!  DETADY, DXSIDX and DXSIDY, and we know that the jacobian gives
!  us the following relation between derivatives with respect to
!  XSI and ETA, and derivatives with respect to X and Y:
!
!    d F(X,Y)/dX     (d XSI/dX  d ETA/dX )   ( d F(XSI, ETA)/d XSI )
!    d F(X,Y)/dY  =  (d XSI/dY  d ETA/dY ) * ( d F(XSI, ETA)/d ETA )
!
!  Here is a graph of the (XSI, ETA) reference triangle we will
!  use.
!
!        ^
!        |
!      1 +        2
!        |       /|
!  ETA   |      4 5
!        |     /  |
!      0 +    1-6-3
!        |
!        +----+---+--->
!             0   1
!
!            XSI
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DET    Output, real ( kind = 8 ) DET, the determinant of the jacobian
!         of the transformation between the reference and isoparametric
!         elements.
!
!  DETADX,
!  DETADY Output, real ( kind = 8 ) DETADX, DETADY, the partial
!         derivative d ETA/d X and d ETA/d Y at (XSI,ETA).
!
!  DXSIDX,
!  DXSIDY Output, real ( kind = 8 ) DXSIDX, DXSIDY, the partial
!         derivative d XSI/d X and d XSI/d Y at (XSI,ETA).
!
!  ETA    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!  IELEM  Input, integer ( kind = 4 ) IELEM, the number of the isoparametric
!         element we are examining.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,NELEM), contains the numbers
!         of the nodes that make up each element.  Element number
!         I is associated with nodes NODE(1,I) through NODE(6,I).
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  XC     Input, real ( kind = 8 ) XC(NP), the X coordinates of the
!         nodes.
!
!  XSI    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!  YC     Input, real ( kind = 8 ) YC(NP), the Y coordinates of the
!         nodes.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) det
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dxdeta
  real ( kind = 8 ) dxdxsi
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) dydeta
  real ( kind = 8 ) dydxsi
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ) eta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) x
  real ( kind = 8 ) xn(6)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) y
  real ( kind = 8 ) yn(6)
  real ( kind = 8 ) yc(np)
!
!  Pick off the X, Y coordinates of the nodes and store them
!  in two short lists.
!
  do i = 1,6
    xn(i) = xc(node(i,ielem))
    yn(i) = yc(node(i,ielem))
  end do
!
!  Set the coefficients in the transformation
!
!    (XSI,ETA) --> (X,Y).
!
!  The mapping has the form:
!
!    X(XSI,ETA) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
!               + D1 * XSI    + E1 * ETA     + F1
!
!    Y(XSI,ETA) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
!               + D2 * XSI    + E2 * ETA     + F2
!
  a1 =  2.0D+00 *xn(1)+2.0D+00 *xn(3)-4.0D+00 *xn(6)
  b1 = -4.0D+00 *xn(3)-4.0D+00 *xn(4)+4.0D+00 *xn(5)+4.0D+00 *xn(6)
  c1 =  2.0D+00 *xn(2)+2.0D+00 *xn(3)-4.0D+00 *xn(5)
  d1 = -3.0D+00 *xn(1)      -xn(3)+4.0D+00 *xn(6)
  e1 =       -xn(2)      +xn(3)+4.0D+00 *xn(4)-4.0D+00 *xn(6)
  f1 =        xn(1)

  a2 =  2.0D+00 *yn(1)+2.0D+00 *yn(3)-4.0D+00 *yn(6)
  b2 = -4.0D+00 *yn(3)-4.0D+00 *yn(4)+4.0D+00 *yn(5)+4.0D+00 *yn(6)
  c2 =  2.0D+00 *yn(2)+2.0D+00 *yn(3)-4.0D+00 *yn(5)
  d2 = -3.0D+00 *yn(1)      -yn(3)+4.0D+00 *yn(6)
  e2 =       -yn(2)      +yn(3)+4.0D+00 *yn(4)-4.0D+00 *yn(6)
  f2 =        yn(1)
!
!  Compute the partial derivatives at the point (XSI,ETA).
!  This is the jacobian matrix
!
!    J: (XSI,ETA) --> (X,Y).
!
  dxdxsi =  2.0D+00 *a1*xsi +       b1*eta + d1
  dxdeta =        b1*xsi + 2.0D+00 *c1*eta + e1

  dydxsi =  2.0D+00 *a2*xsi +       b2*eta + d2
  dydeta =        b2*xsi + 2.0D+00 *c2*eta + e2
!
!  Compute the determinant of the jacobian matrix:
!
!    J: (XSI,ETA) --> (X,Y)
!
  det = dxdxsi*dydeta-dxdeta*dydxsi
!
!  Watch out for a zero determinant.
!
  if ( det == 0.0D+00 ) then
    write ( *, * ) ' '
    write ( *, * ) 'Trans - Fatal error!'
    write ( *, * ) '  The jacobian J: (XSI,ETA) --> (X,Y) is singular!'
    write ( *, * ) '  This occurred for element number ',ielem
    write ( *, * ) '  Local coordinates XSI,ETA = ',xsi,eta
    x = a1*xsi**2+b1*xsi*eta+c1*eta**2+d1*xsi+e1*eta+f1
    y = a2*xsi**2+b2*xsi*eta+c2*eta**2+d2*xsi+e2*eta+f2
    write ( *, * ) '  Global coordinates X,Y = ',x,y
    write ( *, * ) ' '
    write ( *, * ) '  The X, Y nodes were:'
    write ( *, * ) ' '
    do i = 1,6
      write ( *, * ) xn(i),yn(i)
    end do

    stop
  end if
!
!  Compute
!
!    d ETA/d X, d ETA/d Y, d XSI/d X, d XSI/d Y
!
!  by inverting the jacobian matrix
!
!    J: (XSI,ETA) --> (X,Y)
!
!  to get the jacobian matrix
!
!    J: (X,Y) --> (XSI,ETA).
!
!  This uses the simple fact that the inverse of
!
!    (a b)
!    (c d)
!
!  is
!
!    1/(ad-bc) * ( d -b)
!                (-c  a)
!
  dxsidx =  dydeta/det
  dxsidy = -dxdeta/det

  detadx = -dydxsi/det
  detady =  dxdxsi/det

  return
end
subroutine uvpfl(detadx,detady,dpdx,dpdy,dudx,dudy,dvdx,dvdy, &
  dxsidx,dxsidy,eta,gfl,ielem,indx,isotri,nelem,neqnfl,node, &
  np,p,u,v,xc,xq,xsi,yc,yq)
!
!*****************************************************************************80
!
!! UVPFL evaluates the velocities and pressure, and their X and Y
!  derivatives at an arbitrary point in a given element, given
!  the finite element coefficients that represent this data.
!
!  If the element is not isoparametric, then UVPFL requires the
!  physical X and Y coordinates of the point.
!
!  If the element is isoparametric, UVPFL requires the XSI, ETA
!  coordinates of the point. 
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DETADX,
!  DETADY Input, real ( kind = 8 ) DETADX, DETADY, the partial derivative
!         d ETA/d X and d ETA/d Y at (XSI,ETA).
!
!  DPDX,
!  DPDY   Output, real ( kind = 8 ) DPDX, DPDY, the partial derivatives
!         of P with respect to X and Y.
!
!  DUDX,
!  DUDY   Output, real ( kind = 8 ) DUDX, DUDY, the partial derivatives
!         of U with respect to X and Y.
!
!  DVDX,
!  DVDY   Output, real ( kind = 8 ) DVDX, DVDY, the partial derivatives
!         of V with respect to X and Y.
!
!  DXSIDX,
!  DXSIDY Input, real ( kind = 8 ) DXSIDX, DXSIDY, the partial derivative
!         d XSI/d X and d XSI/d Y at (XSI,ETA).
!
!  ETA    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point,
!         needed only if the element is isoparametric.
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!
!         GFL is the computed solution vector, in which are stored
!         pressures and velocities.
!
!  IELEM  Input, integer ( kind = 4 ) IELEM, the element in which the point lies
!         at which the quantities are desired.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  ISOTRI Input, integer ( kind = 4 ) ISOTRI(NELEM).
!
!         0, the element is NOT isoparametric, and the nodes never move.
!         That means that the quadrature points are only computed once.
!
!         1, the element is NOT isoparametric, but the nodes may move.
!         Quadrature point locations must be updated on each step.
!         This could occur for elements above, but not touching, the bump.
!
!         2, the element is isoparametric.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of equations in the full system.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,MAXELM).
!
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  P      Output, real ( kind = 8 ) P, the pressure.
!
!  U      Output, real ( kind = 8 ) U, the horizontal velocity.
!
!  V      Output, real ( kind = 8 ) V, the vertical velocity.
!
!  XC     Input, real ( kind = 8 ) XC(NP), the X coordinates of the nodes.
!
!  XQ     Input, real ( kind = 8 ) XQ, the X coordinate of the point.
!
!  XSI    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point,
!         needed only if the element is isoparametric.
!
!  YC     Input, real ( kind = 8 ) YC(NP), the Y coordinates of the nodes.
!
!  YQ     Input, real ( kind = 8 ) YQ, the Y coordinate of the point.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
!
  real ( kind = 8 ) coef
  real ( kind = 8 ) detadx
  real ( kind = 8 ) detady
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) dxsidx
  real ( kind = 8 ) dxsidy
  real ( kind = 8 ) eta
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) iun
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p
  real ( kind = 8 ) q
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xq
  real ( kind = 8 ) xsi
  real ( kind = 8 ) yc(np)
  real ( kind = 8 ) yq
!
  p = 0.0D+00
  u = 0.0D+00
  v = 0.0D+00
  dpdx = 0.0D+00
  dpdy = 0.0D+00
  dudx = 0.0D+00
  dudy = 0.0D+00
  dvdx = 0.0D+00
  dvdy = 0.0D+00

  do iq = 1,6
!
!  Evaluate the basis functions W and Q, and their derivatives
!  DQDX, DQDY, DWDX, DWDY at XQ, YQ.
!
    if ( isotri(ielem) == 0.or.isotri(ielem).eq.1) then

      call qbf(ielem,iq,w,dwdx,dwdy,nelem,node,np,xc,xq,yc,yq)

      call bsp(q,dqdx,dqdy,ielem,iq,nelem,node,np,xc,xq,yc,yq)

    else

      call refqbf(w,dwdx,dwdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)

      call refbsp(q,dqdx,dqdy,detadx,detady,dxsidx,dxsidy,eta,iq,xsi)

    end if
!
!  Compute the coefficients at the node at XP, YP.
!
    ip = node(iq,ielem)

    iun = indx(1,ip)
    coef = gfl(iun)
    u = u+coef*w
    dudx = dudx+coef*dwdx
    dudy = dudy+coef*dwdy

    iun = indx(2,ip)
    coef = gfl(iun)
    v = v+coef*w
    dvdx = dvdx+coef*dwdx
    dvdy = dvdy+coef*dwdy

    iun = indx(3,ip)
    if ( 0 < iun ) then
      coef = gfl(iun)
      p = p+coef*q
      dpdx = dpdx+coef*dqdx
      dpdy = dpdy+coef*dqdy
    end if

  end do

  return
end
subroutine uvpnrm(gfl,indx,neqnfl,np,pnorm,uvnorm)
!
!*****************************************************************************80
!
!! UVPNRM returns the "norm" of the solution.  Here, the norm of
!  a solution GFL = (U,V,P) is defined as two numbers, the maximum
!  velocity magnitude at a node, and the maximum pressure at a node.
!
! 
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL).
!
!         GFL is the current solution vector, in which are stored
!         the finite element coefficients that define the velocity
!         and pressure functions, U, V and P.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of finite element equations used
!         to define the horizontal and vertical velocities and the
!         pressure.
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  PNORM  Output, real ( kind = 8 ) PNORM, the maximum absolute value
!         pressure coefficient.
!
!  UVNORM Output, real ( kind = 8 ) UVNORM, the maximum velocity magnitude.
!
  implicit none

  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np

  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indx(3,np)
  real ( kind = 8 ) p
  real ( kind = 8 ) pnorm
  real ( kind = 8 ) u
  real ( kind = 8 ) uvnorm
  real ( kind = 8 ) v

  uvnorm = 0.0D+00
  pnorm = 0.0D+00

  do i = 1,np

    u = gfl(indx(1,i))
    v = gfl(indx(2,i))
    uvnorm = max(uvnorm,sqrt(u**2+v**2))

    if ( 0 < indx(3,i) ) then
      p = gfl(indx(3,i))
      pnorm = max(pnorm,abs(p))
    end if

  end do

  return
end
subroutine uvpqfl(dpdx,dpdy,dudx,dudy,dvdx,dvdy,gfl,ielem,indx, &
  iquad,nelem,neqnfl,node,np,p,phifl,u,v)
!
!*****************************************************************************80
!
!! UVPQFL evaluates the velocities and pressure, and their X and Y
!  derivatives, at a quadrature point in a given element, given
!  the finite element coefficients, and the value of the basis functions
!  and their X and Y derivatives at each quadrature point.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DPDX,
!  DPDY   Output, real ( kind = 8 ) DPDX, DPDY, the derivatives of the
!         pressure function with respect to X and Y.
!
!  DUDX,
!  DUDY   Output, real ( kind = 8 ) DUDX, DUDY, the derivatives of the
!         horizontal velocity function with respect to X and Y.
!
!  DVDX,
!  DVDY   Output, real ( kind = 8 ) DVDX, DVDY, the derivatives of the
!         vertical velocity function with respect to X and Y.
!
!  GFL    Input, real ( kind = 8 ) GFL(NEQNFL), the current solution
!         estimate for the full problem.
!
!  IELEM  Input, integer ( kind = 4 ) IELEM, the element in which the quadrature
!         point lies.
!
!  INDX   Input, integer ( kind = 4 ) INDX(3,NP). 
!
!         INDX(I,J) contains, for each node J, the index of U, V and P at
!         that node, or 0 or a negative value.
!
!         If K = INDX(I,J) is positive, then the value of the degree
!         of freedom is stored in the solution vector entry GFL(K),
!         and an equation will be generated to determine its value.
!
!         If INDX(I,J) is not positive, then no equation is
!         generated to determine for variable I at node J, either because
!         the variable is specified in some other way, or because
!         (in the case of pressure), there is no coefficient associated
!         with that node.
!
!  IQUAD  Input, integer ( kind = 4 ) IQUAD, the local index of the quadrature point.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NEQNFL Input, integer ( kind = 4 ) NEQNFL, the number of equations in the full system.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,MAXELM).
!
!         NODE(I,J) contains, for an element J, the global node index of
!         the element node whose local number is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!               2
!              /|
!             4 5
!            /  |
!           1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes used to define the finite
!         element mesh.  NP = (2*NX-1)*(2*NY-1).
!
!  P      Output, real ( kind = 8 ) P, the value of the pressure.
!
!  PHIFL  Input, real ( kind = 8 ) PHIFL(3,6,10,NELEM). 
!
!         PHIFL contains the value of a finite element basis function, its
!         derivative, or other information, evaluated at the quadrature
!         points.
!
!         The meaning of the entry PHIFL(I,J,K,L) is as follows. 
!         For the quadrature point I, and basis function J, in element L,
!         PHIFL(I,J,K,L) represents the value of:
!
!           K =  1, W, the finite element basis function for velocities;
!           K =  2, dWdX, the X derivative of W;
!           K =  3, dWdY, the Y derivative of W;
!           K =  4, Q, the finite element basis function for pressures;
!           K =  5, dQdX, the X derivative of Q;
!           K =  6, dQdY, the Y derivative of Q;
!           K =  7, dXsidX, the X derivative of the mapping (X,Y)->XSI;
!           K =  8, dXsidY, the Y derivative of the mapping (X,Y)->XSI;
!           K =  9, dEtadX, the X derivative of the mapping (X,Y)->ETA;
!           K = 10, dEtadY, the Y derivative of the mapping (X,Y)->ETA;
!
!         In particular, PHIFL(I,J,K,L) is the value of the quadratic
!         basis function W associated with local node J in element L,
!         evaluated at quadrature point I.
!
!         Note that PHIFL(I,J,K,L) = 0 whenever J=4, 5, or 6 and K=4, 5, or 6,
!         since there are only three linear basis functions.
!
!  U      Output, real ( kind = 8 ) U, the value of the horizontal
!         velocity.
!
!  V      Output, real ( kind = 8 ) V, the value of the vertical
!         velocity.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np

  real ( kind = 8 ) coef
  real ( kind = 8 ) dpdx
  real ( kind = 8 ) dpdy
  real ( kind = 8 ) dqdx
  real ( kind = 8 ) dqdy
  real ( kind = 8 ) dudx
  real ( kind = 8 ) dudy
  real ( kind = 8 ) dvdx
  real ( kind = 8 ) dvdy
  real ( kind = 8 ) dwdx
  real ( kind = 8 ) dwdy
  real ( kind = 8 ) gfl(neqnfl)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iq
  integer ( kind = 4 ) iquad
  integer ( kind = 4 ) iun
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p
  real ( kind = 8 ) phifl(3,6,10,nelem)
  real ( kind = 8 ) q
  real ( kind = 8 ) u
  real ( kind = 8 ) v
  real ( kind = 8 ) w
!
!  Start all the functions at zero.
!
  p = 0.0D+00
  u = 0.0D+00
  v = 0.0D+00
  dpdx = 0.0D+00
  dpdy = 0.0D+00
  dudx = 0.0D+00
  dudy = 0.0D+00
  dvdx = 0.0D+00
  dvdy = 0.0D+00
!
!  Now each of these functions is represented as the sum of
!  coefficients times basis functions.  In this particular
!  element, at this particular quadrature point, we know that
!  exactly 6 basis functions are nonzero.  So if
!  we simply look up the values of the basis functions (and
!  their X and Y derivatives), and multiply by the appropriate
!  coefficients, we can evaluate the functions.
!
!  W, DWDX and DWDY represent the value of a quadratic basis
!  function and its X and Y derivative.
!
!  Q, DQDX and DQDY represent the value of a linear basis
!  function and its X and Y derivatives.
!
  do iq = 1,6

    w = phifl(iquad,iq,1,ielem)
    dwdx = phifl(iquad,iq,2,ielem)
    dwdy = phifl(iquad,iq,3,ielem)

    q = phifl(iquad,iq,4,ielem)
    dqdx = phifl(iquad,iq,5,ielem)
    dqdy = phifl(iquad,iq,6,ielem)
!
!  Now that we have the basis function values, we need to look
!  up the coefficient COEF that multiplies the basis function.
!
    ip = node(iq,ielem)

    iun = indx(1,ip)
    coef = gfl(iun)
    u = u+coef*w
    dudx = dudx+coef*dwdx
    dudy = dudy+coef*dwdy

    iun = indx(2,ip)
    coef = gfl(iun)
    v = v+coef*w
    dvdx = dvdx+coef*dwdx
    dvdy = dvdy+coef*dwdy

    iun = indx(3,ip)
    if ( 0 < iun ) then
      coef = gfl(iun)
      p = p+coef*q
      dpdx = dpdx+coef*dqdx
      dpdy = dpdy+coef*dqdy
    end if

  end do

  return
end
subroutine uvpqrb(dprbdx,dprbdy,durbdx,durbdy,dvrbdx,dvrbdy,grb, &
  ielem,iquad,maxcofrb,maxelm,ncofrb,phirb,prb,urb,vrb)
!
!*****************************************************************************80
!
!! UVPQRB evaluates reduced basis state variables at a quadrature point.
!
!
!  Discussion:
!
!    The routine is given:
!
!      GRB, the reduced basis coefficients,
!      PHIRB, the reduced basis functions and derivatives evaluated
!        at the quadrature points,
!      IELEM, a particular element, and
!      IQUAD, a particular quadrature point in that element.
!
!    and computes:
!
!      URB, VRB, PRB, the state variables evaluated at that quadrature point,
!      the X and Y derivatives of these quantities.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    09 July 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  DPRBDX,
!  DPRBDY Output, real ( kind = 8 ) DPRBDX, DPRBDY, the derivatives of the
!         pressure function with respect to X and Y.
!
!  DURBDX,
!  DURBDY Output, real ( kind = 8 ) DURBDX, DURBDY, the derivatives of the
!         horizontal velocity function with respect to X and Y.
!
!  DVRBDX,
!  DVRBDY Output, real ( kind = 8 ) DVRBDX, DVRBDY, the derivatives of the
!         vertical velocity function with respect to X and Y.
!
!  GRB    Input, real ( kind = 8 ) GRB(NCOFRB), the current solution
!         estimate for the reduced problem.
!
!  IELEM  Input, integer ( kind = 4 ) IELEM, the element in which the quadrature
!         point lies.
!
!  IQUAD  Input, integer ( kind = 4 ) IQUAD, the local index of the quadrature point.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NCOFRB Input, integer ( kind = 4 ) NCOFRB, the number of equations in the reduced
!         system.
!
!  PRB    Output, real ( kind = 8 ) PRB, the value of the pressure.
!
!  PHIRB  Input, real ( kind = 8 ) PHIRB(3,NCOFRB,15,NELEM).
!         PHIRB contains the values of a finite element basis function
!         or its X or Y derivative, in a given element, at a given
!         quadrature point, for a particular reduced basis function.
!
!         For PHIRB(I,J,K,L), index J refers to the reduced basis
!         basis functions, for J = 0 to NCOFRB.
!
!         The meaning of the K index of PHIRB(I,J,K,L) is as follows:
!
!           For the quadrature point I, and reduced basis function J,
!           in element L, PHIRB(I,J,K,L) represents the value of:
!
!             K = 1, WUrb, the finite element U velocity basis function;
!             K = 2, dWUrbdX, the X derivative of WUrb;
!             K = 3, dWUrbdY, the Y derivative of WUrb;
!             K = 4, WVrb, the finite element V velocity basis function;
!             K = 5, dWVrbdX, the X derivative of WVrb;
!             K = 6, dWVrbdY, the Y derivative of WVrb;
!             K = 7, Q, the finite element pressure basis function.
!             K = 8, dQrbdX, the X derivative of Qrb;
!             K = 9, dQrbdY, the Y derivative of Qrb.
!             K = 10, WU0rb, same as WUrb, with zero BC.
!             K = 11, dWU0rbdX, same as dWUrbdX, with zero BC.
!             K = 12, dWU0rbdY, same as dWUrbdY, with zero BC.
!             K = 13, WV0rb, same as WVrb, with zero BC.
!             K = 14, dWV0rbdX, same as dWVrbdX, with zero BC.
!             K = 15, dWV0rbdY, same as dWVrbdY, with zero BC.
!
!  URB    Output, real ( kind = 8 ) URB, the value of the horizontal
!         velocity.
!
!  VRB    Output, real ( kind = 8 ) VRB, the value of the vertical
!         velocity.
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxelm
  integer ( kind = 4 ) ncofrb
!
  real ( kind = 8 ) coef
  real ( kind = 8 ) dprbdx
  real ( kind = 8 ) dprbdy
  real ( kind = 8 ) dqrbdx
  real ( kind = 8 ) dqrbdy
  real ( kind = 8 ) durbdx
  real ( kind = 8 ) durbdy
  real ( kind = 8 ) dvrbdx
  real ( kind = 8 ) dvrbdy
  real ( kind = 8 ) dwurbdx
  real ( kind = 8 ) dwurbdy
  real ( kind = 8 ) dwvrbdx
  real ( kind = 8 ) dwvrbdy
  real ( kind = 8 ) grb(ncofrb)
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) ieqnrb
  integer ( kind = 4 ) iquad
  real ( kind = 8 ) prb
  real ( kind = 8 ) phirb(3,maxcofrb,15,maxelm)
  real ( kind = 8 ) qrb
  real ( kind = 8 ) urb
  real ( kind = 8 ) vrb
  real ( kind = 8 ) wurb
  real ( kind = 8 ) wvrb
!
!  Start all the function values at 0.
!
  urb = 0.0D+00
  durbdx = 0.0D+00
  durbdy = 0.0D+00

  vrb = 0.0D+00
  dvrbdx = 0.0D+00
  dvrbdy = 0.0D+00

  prb = 0.0D+00
  dprbdx = 0.0D+00
  dprbdy = 0.0D+00
!
!  Now each of these functions is represented as the sum of
!  coefficients times the reduced basis vectors.  So if
!  we simply look up the values of the reduced basis functions (and
!  their X and Y derivatives), and multiply by the appropriate
!  coefficients, we can evaluate the functions.
!
!  WURB, DWURBDX and DWURBDY represent the value of a quadratic basis
!  function and its X and Y derivative.
!
!  WVRB, DWVRBDX and DWVRBDY represent the value of a quadratic basis
!  function and its X and Y derivative.
!
!  QRB, DQRBDX and DQRBDY represent the value of a linear basis
!  function and its X and Y derivatives.
!
!  See routine SETPRB, where these values are loaded into PHIRB.
!
  do ieqnrb = 1,ncofrb

    wurb    = phirb(iquad,ieqnrb,1,ielem)
    dwurbdx = phirb(iquad,ieqnrb,2,ielem)
    dwurbdy = phirb(iquad,ieqnrb,3,ielem)

    wvrb    = phirb(iquad,ieqnrb,4,ielem)
    dwvrbdx = phirb(iquad,ieqnrb,5,ielem)
    dwvrbdy = phirb(iquad,ieqnrb,6,ielem)

    qrb    = phirb(iquad,ieqnrb,7,ielem)
    dqrbdx = phirb(iquad,ieqnrb,8,ielem)
    dqrbdy = phirb(iquad,ieqnrb,9,ielem)

    coef = grb(ieqnrb)

    urb    = urb+coef*wurb
    durbdx = durbdx+coef*dwurbdx
    durbdy = durbdy+coef*dwurbdy

    vrb    = vrb+coef*wvrb
    dvrbdx = dvrbdx+coef*dwvrbdx
    dvrbdy = dvrbdy+coef*dwvrbdy

    prb    = prb+coef*qrb
    dprbdx = dprbdx+coef*dqrbdx
    dprbdy = dprbdy+coef*dqrbdy

  end do

  return
end
subroutine wrdis(disfil,eqn,indx,isotri,maxcofrb,maxnfl,ncofrb,nelem, &
  neqnfl,node,np,npar,npe,nprof,nx,ny,p,par,senfl,u,v,xc,xprof,yc)
!
!*****************************************************************************80
!
!! WRDIS writes information to a file which can be used to create
!  graphics images. 
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
  implicit none
!
  integer ( kind = 4 ) maxcofrb
  integer ( kind = 4 ) maxnfl
  integer ( kind = 4 ) ncofrb
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) neqnfl
  integer ( kind = 4 ) np
  integer ( kind = 4 ) npar
  integer ( kind = 4 ) ny
!
  character ( len = 2 ) ctemp
  character ( len = 30 ) disfil
  character ( len = 2 ) eqn(neqnfl)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) icheck
  integer ( kind = 4 ) igunit
  integer ( kind = 4 ) ihor
  integer ( kind = 4 ) indx(3,np)
  integer ( kind = 4 ) iprs
  integer ( kind = 4 ) isen
  integer ( kind = 4 ) isotri(nelem)
  integer ( kind = 4 ) iver
  integer ( kind = 4 ) j
  integer ( kind = 4 ) node(6,nelem)
  integer ( kind = 4 ) npe
  integer ( kind = 4 ) nprof(2*ny-1)
  integer ( kind = 4 ) nsen
  integer ( kind = 4 ) nx
  real ( kind = 8 ) p(np)
  real ( kind = 8 ) par(npar)
  real rtemp
  real ( kind = 8 ) senfl(maxnfl,maxcofrb)
  real ( kind = 8 ) u(np)
  real ( kind = 8 ) v(np)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xprof
  real ( kind = 8 ) yc(np)
!
!  Delete any old copy of the file.
!
  igunit = 11

  open(unit = igunit,file=disfil,status='unknown', &
    form = 'formatted',access='sequential')

  nsen = ncofrb
!
!  Number of elements, nodes, parameters,
!  elements in the X direction, elements in the Y direction.
!
  write(igunit,*)nelem
  write(igunit,*)np
  write(igunit,*)npar
  write(igunit,*)npe
  write(igunit,*)nsen
  write(igunit,*)nx
  write(igunit,*)ny
!
!  Pressures, P.
!
  do i = 1,np
    write(igunit,*)p(i)
  end do
!
!  Horizontal velocities, U.
!
  do i = 1,np
    write(igunit,*)u(i)
  end do
!
!  Vertical velocities, V
!
  do i = 1,np
    write(igunit,*)v(i)
  end do
!
!  Indicator of element type (isoparametric or not).
!
  do i = 1,nelem
    write(igunit,*)isotri(i)
  end do
!
!  Nodes that make up each element.
!
  do i = 1,npe
    do j = 1,nelem
      write(igunit,*)node(i,j)
    end do
  end do
!
!  Indices of the nodes along the profile line.
!
  do i = 1,2*ny-1
    write(igunit,*)nprof(i)
  end do
!
!  Parameters.
!
  do i = 1,npar
    write(igunit,*)par(i)
  end do
!
!  Pressure sensitivities, dP/dpar
!
  do isen = 1,nsen

    do i = 1,np
      iprs = indx(3,i)
      rtemp = 0.0D+00
      if ( 0 < iprs ) then
        rtemp = sngl(senfl(iprs,isen))
      end if
      write(igunit,*)rtemp
    end do
!
!  Horizontal velocity sensitivities, dU/dpar
!
    do i = 1,np
      ihor = indx(1,i)
      write(igunit,*)senfl(ihor,isen)
    end do
!
!  Vertical velocity sensitivities, dV/dpar
!
    do i = 1,np
      iver = indx(2,i)
      write(igunit,*)senfl(iver,isen)
    end do

  end do
!
!  X coordinates of nodes.
!
  do i = 1,np
    write(igunit,*)xc(i)
  end do
!
!  X coordinate of profile line.
!
  write(igunit,*)xprof
!
!  Y coordinates of nodes.
!
  do i = 1,np
    write(igunit,*)yc(i)
  end do
!
!  Nodal equation types.
!
  do i = 1,np
    ihor = indx(1,i)
    iver = indx(2,i)
    iprs = indx(3,i)
    if ( iprs <= 0) then
      ctemp = '  '
    else
      ctemp = eqn(iprs)
    end if
    write(igunit,'(3a2)')eqn(ihor),eqn(iver),ctemp
  end do
!
!  Write a check at the end.
!
  icheck = 1953
  write(igunit,*)icheck

  close(unit = igunit)

  return
end
subroutine wrtec(nelem,node,np,p,tecfil,u,v,xc,yc)
!
!*****************************************************************************80
!
!! WRTEC writes data suitable for plotting by TECPLOT.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 February 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,MAXELM) or NODE(6,NELEM).
!
!         NODE(I,J) contains, for an element J, the global index of
!         the node whose local number in J is I.
!
!         The local ordering of the nodes is suggested by this diagram:
!
!           Global nodes   Elements      NODE
!                                                          1  2  3  4  5  6
!           74  84  94     3-6-1   2     Left element =  (94,72,74,83,73,84)
!                          |  /   /|
!           73  83  93     5 4   4 5     Right element = (72,94,92,83,93,82)
!                          |/   /  |
!           72  82  92     2   1-6-3
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  P      Input, real P(NP), the pressure.
!
!  U      Input, real U(NP), the horizontal velocity.
!
!  V      Input, real V(NP), the vertical velocity.
!
!  XC     Input, real XC(NP), the X coordinates of the nodes.
!
!  YC     Input, real YC(NP), the Y coordinates of the nodes.
!
  implicit none
!
  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) p(np)
  character ( len = 30 ) tecfil
  real ( kind = 8 ) u(np)
  real ( kind = 8 ) v(np)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) yc(np)
!
!  Delete any old copy of the TECPLOT data file.
!
  open(unit = 10,file=tecfil,status='old',err=10)
  write ( *, * ) ' '
  write ( *, * ) 'WRTEC - Note:'
  write ( *, * ) '  Deleting an old copy of the TECPLOT data file.'
  close(unit = 10,status='delete')

10    continue

  open(unit = 10,file=tecfil,status='unknown')

  write(10,*)'Title = ','"ARBY data"'
  write(10,*)'Variables = "X","Y","P","U","V"'
  write(10,*)'Zone N = ',np,', E=',4*nelem,', F=FEPOINT, ET=TRIANGLE'
!
!  Write out the data at each node.
!
  do i = 1,np
    write(10,'(5g15.6)')xc(i),yc(i),p(i),u(i),v(i)
  end do
!
!  Write out the data that defines the elements.
!  Each 6 node quadratic element must be described as 4 linear elements.
!
  do i = 1,nelem
    write(10,'(3i6)')node(1,i),node(4,i),node(6,i)
    write(10,'(3i6)')node(2,i),node(5,i),node(4,i)
    write(10,'(3i6)')node(3,i),node(6,i),node(5,i)
    write(10,'(3i6)')node(4,i),node(5,i),node(6,i)
  end do

  close(unit = 10)

  return
end
subroutine xofxsi(eta,ielem,nelem,node,np,x,xc,xsi,y,yc)
!
!*****************************************************************************80
!
!! XOFXSI is given the XSI, ETA coordinates of a point in an
!  isoparametric element and determines its X, Y coordinates.
!
!  Here is a graph of the (XSI, ETA) reference triangle we will use.
!
!        ^
!        |
!      1 +        2
!        |       /|
!  ETA   |      4 5
!        |     /  |
!      0 +    1-6-3
!        |
!        +----+---+--->
!             0   1
!
!            XSI
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 October 1996
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!  ETA    Input, real ( kind = 8 ) ETA, the ETA coordinate of the point.
!
!  IELEM  Input, integer ( kind = 4 ) IELEM, the number of the isoparametric
!         element we are examining.
!
!  NELEM  Input, integer ( kind = 4 ) NELEM, the number of elements.
!
!  NODE   Input, integer ( kind = 4 ) NODE(6,nelem), contains the numbers
!         of the nodes that make up each element.  Element number
!         I is associated with nodes NODE(1,I) through NODE(6,I).
!
!  NP     Input, integer ( kind = 4 ) NP, the number of nodes.
!
!  X      Output, real ( kind = 8 ) X, the X coordinate of the point.
!
!  XC     Input, real ( kind = 8 ) XC(NP), the X coordinates of the
!         nodes.
!
!  XSI    Input, real ( kind = 8 ) XSI, the XSI coordinate of the point.
!
!  Y      Output, real ( kind = 8 ) Y, the Y coordinate of the point.
!
!  YC     Input, real ( kind = 8 ) YC(NP), the Y coordinates of the
!         nodes.
!
  implicit none

  integer ( kind = 4 ) nelem
  integer ( kind = 4 ) np

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  real ( kind = 8 ) d1
  real ( kind = 8 ) d2
  real ( kind = 8 ) e1
  real ( kind = 8 ) e2
  real ( kind = 8 ) eta
  real ( kind = 8 ) f1
  real ( kind = 8 ) f2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ielem
  integer ( kind = 4 ) node(6,nelem)
  real ( kind = 8 ) x
  real ( kind = 8 ) xn(6)
  real ( kind = 8 ) xc(np)
  real ( kind = 8 ) xsi
  real ( kind = 8 ) y
  real ( kind = 8 ) yn(6)
  real ( kind = 8 ) yc(np)
!
!  Pick off the X, Y coordinates of the nodes and store them
!  in two short lists.
!
  do i = 1,6
    xn(i) = xc(node(i,ielem))
    yn(i) = yc(node(i,ielem))
  end do
!
!  Set the coefficients in the transformation
!
!    (XSI,ETA) --> (X,Y).
!
!  The mapping has the form:
!
!    X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
!               + D1 * XSI    + E1 * ETA     + F1
!
!    Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
!               + D2 * XSI    + E2 * ETA     + F2
!
  a1 =  2.0D+00 *xn(1)+2.0D+00 *xn(3)-4.0D+00 *xn(6)
  b1 = -4.0D+00 *xn(3)-4.0D+00 *xn(4)+4.0D+00 *xn(5)+4.0D+00 *xn(6)
  c1 =  2.0D+00 *xn(2)+2.0D+00 *xn(3)-4.0D+00 *xn(5)
  d1 = -3.0D+00 *xn(1)      -xn(3)+4.0D+00 *xn(6)
  e1 =       -xn(2)      +xn(3)+4.0D+00 *xn(4)-4.0D+00 *xn(6)
  f1 =        xn(1)

  a2 =  2.0D+00 *yn(1)+2.0D+00 *yn(3)-4.0D+00 *yn(6)
  b2 = -4.0D+00 *yn(3)-4.0D+00 *yn(4)+4.0D+00 *yn(5)+4.0D+00 *yn(6)
  c2 =  2.0D+00 *yn(2)+2.0D+00 *yn(3)-4.0D+00 *yn(5)
  d2 = -3.0D+00 *yn(1)      -yn(3)+4.0D+00 *yn(6)
  e2 =       -yn(2)      +yn(3)+4.0D+00 *yn(4)-4.0D+00 *yn(6)
  f2 =        yn(1)

  x = a1*xsi**2 + b1*xsi*eta + c1*eta**2 + d1*xsi + e1*eta + f1

  y = a2*xsi**2 + b2*xsi*eta + c2*eta**2 + d2*xsi + e2*eta + f2

  return
end
function dasum(n,dx,incx)
!
!*****************************************************************************80
!
!! DASUM takes the sum of the absolute values of the entries of
!  a vector.
!
  real ( kind = 8 ) dasum
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
!
  if ( n <= 0) then
    dasum = 0.0d0
    return
  end if

  if ( incx <= 0) then

    dasum = 0.0d0
    return

  else if ( incx /= 1) then

    dtemp = 0.0d0
    do i = 1,n*incx,incx
      dtemp = dtemp+abs(dx(i))
    end do

  else

    m = mod(n,6)

    dtemp = 0.0d0
    do i = 1,m
      dtemp = dtemp+abs(dx(i))
    end do

    do i = m+1,n,6
      dtemp = dtemp+abs(dx(i))+abs(dx(i+1))+abs(dx(i+2)) &
        +abs(dx(i+3))+abs(dx(i+4))+abs(dx(i+5))
    end do

  end if

  dasum = dtemp

  return
end
subroutine daxpy(n,da,dx,incx,dy,incy)
!
!***********************************************************************
!
!! DAXPY adds a multiple of one vector to another.
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) da
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
!
  if ( n <= 0)return

  if ( da == 0.0D+00 )return

  if ( incx == 1.and.incy == 1)go to 20

  ix = 1
  iy = 1
  if ( incx < 0)ix = (-n+1)*incx + 1
  if ( incy < 0)iy = (-n+1)*incy + 1

  do i = 1,n
    dy(iy) = dy(iy) + da*dx(ix)
    ix = ix + incx
    iy = iy + incy
  end do

  return

   20 m = mod(n,4)

  do i = 1,m
    dy(i) = dy(i) + da*dx(i)
  end do

  do i = m+1,n,4
    dy(i) = dy(i) + da*dx(i)
    dy(i + 1) = dy(i + 1) + da*dx(i + 1)
    dy(i + 2) = dy(i + 2) + da*dx(i + 2)
    dy(i + 3) = dy(i + 3) + da*dx(i + 3)
  end do

  return
end
subroutine dcopy(n,dx,incx,dy,incy)
!
!***********************************************************************
!
!! DCOPY copies a vector X to a vector Y.
!
!
!  N      Input, integer ( kind = 4 ) N, the number of entries to copy.
!
!  DX     Input, real ( kind = 8 ) DX(*), the vector to be copied.
!
!  INCX   Input, integer ( kind = 4 ) INCX, the increment between successive
!         entries of DX.
!
!  DY     Output, real ( kind = 8 ) DY(*), the vector to be copied.
!
!  INCY   Input, integer ( kind = 4 ) INCY, the increment between successive
!         entries of DY.
!
  implicit none
!
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
!
  if ( n <= 0)return

  if ( incx == 1.and.incy == 1) then

    m = mod(n,7)

    do i = 1,m
      dy(i) = dx(i)
    end do

    do i = m+1,n,7
      dy(i) = dx(i)
      dy(i + 1) = dx(i + 1)
      dy(i + 2) = dx(i + 2)
      dy(i + 3) = dx(i + 3)
      dy(i + 4) = dx(i + 4)
      dy(i + 5) = dx(i + 5)
      dy(i + 6) = dx(i + 6)
    end do

  else

    if ( incx < 0) then
      ix = (-n+1)*incx + 1
    else
      ix = 1
    end if

    if ( incy < 0) then
      iy = (-n+1)*incy + 1
    else
      iy = 1
    end if

    do i = 1,n
      dy(iy) = dx(ix)
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
function ddot(n,dx,incx,dy,incy)
!
!***********************************************************************
!
!! DDOT forms the dot product of two vectors.
!
  real ( kind = 8 ) ddot
  real ( kind = 8 ) dtemp
  real ( kind = 8 ) dx(*)
  real ( kind = 8 ) dy(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
!
  ddot = 0.0d0
  dtemp = 0.0d0

  if ( n <= 0)return

  if ( incx == 1.and.incy == 1)go to 20

  if ( incx < 0) then
    ix = (-n+1)*incx + 1
  else
    ix = 1
  end if

  if ( incy < 0) then
    iy = (-n+1)*incy + 1
  else
    iy = 1
  end if

  do i = 1,n
    dtemp = dtemp + dx(ix)*dy(iy)
    ix = ix + incx
    iy = iy + incy
  end do

  ddot = dtemp
  return
!
   20 m = mod(n,5)

  do i = 1,m
    dtemp = dtemp + dx(i)*dy(i)
  end do

  do i = m+1,n,5
    dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
      dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
  end do

  ddot = dtemp

  return
end
subroutine dgbtf2(m,n,kl,ku,ab,ldab,ipiv,info)
!
!***********************************************************************
!
!! DGBTF2 ???
!
  integer ( kind = 4 )            info, kl, ku, ldab, m, n
!     ..
!     .. Array Arguments ..
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   ab( ldab, * )
!
!  DGBTF2 computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix A.
!
!  N       (input) integer
!          The number of columns of the matrix A.
!
!  KL      (input) integer
!          The number of subdiagonals within the band of A.
!
!  KU      (input) integer
!          The number of superdiagonals within the band of A.
!
!  AB      (input/output) real ( kind = 8 ) array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) integer
!          The leading dimension of the array AB.  2*KL+KU+1 <= LDAB.
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U, because of fill-in resulting from the row
!  interchanges.
!
  real ( kind = 8 )   one, zero
  parameter          ( one = 1.0d+0, zero = 0.0d+0 )
  integer ( kind = 4 )            i, j, jp, ju, km, kv
  integer ( kind = 4 )            idamax
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
  kv = ku + kl
!
!     Test the input parameters.
!
  info = 0
  if ( m < 0 ) then
     info = -1
  else if ( n < 0 ) then
     info = -2
  else if ( kl < 0 ) then
     info = -3
  else if ( ku < 0 ) then
     info = -4
  else if ( ldab < kl+kv+1 ) then
     info = -6
  end if
  if ( info /= 0 ) then
     call xerbla( 'dgbtf2', -info )
     return
  end if

  if ( m == 0 .or. n == 0 ) then
    return
  end if
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
  do j = ku + 2, min( kv, n )
    do i = kv - j + 2, kl
      ab( i, j ) = zero
    end do
  end do
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
  ju = 1
!
  do j = 1, min( m, n )
!
!        Set fill-in elements in column J+KV to zero.
!
     if ( j+kv <= n ) then
        do i = 1, kl
           ab( i, j+kv ) = zero
        end do
     end if
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
     km = min( kl, m-j )
     jp = idamax( km+1, ab( kv+1, j ), 1 )
     ipiv( j ) = jp + j - 1
     if ( ab( kv+jp, j ) /= zero ) then
        ju = max( ju, min( j+ku+jp-1, n ) )
!
!  Apply interchange to columns J to JU.
!
        if ( jp /= 1 ) then
          call dswap( ju-j+1, ab( kv+jp, j ), ldab-1, ab( kv+1, j ), ldab-1 )
        end if

        if ( 0 < km ) then
!
!  Compute multipliers.
!
           call dscal( km, one / ab( kv+1, j ), ab( kv+2, j ), 1 )
!
!  Update trailing submatrix within the band.
!
           if ( j < ju ) then

              call dger( km, ju-j, -one, ab( kv+2, j ), 1, &
                ab( kv, j+1 ), ldab-1, ab( kv+1, j+1 ), ldab-1 )
           end if

        end if
     else
!
!  If pivot is zero, set INFO to the index of the pivot
!  unless a zero pivot has already been found.
!
        if ( info == 0 ) then
           info = j
        end if

     end if

  end do

  return
end
subroutine dgbtrf(m,n,kl,ku,ab,ldab,ipiv,info)
!
!***********************************************************************
!
!! DGBTRF ???
!
  integer ( kind = 4 )            info, kl, ku, ldab, m, n
!     ..
!     .. Array Arguments ..
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   ab( ldab, * )
!
!  DGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix A.
!
!  N       (input) integer
!          The number of columns of the matrix A.
!
!  KL      (input) integer
!          The number of subdiagonals within the band of A.
!
!  KU      (input) integer
!          The number of superdiagonals within the band of A.
!
!  AB      (input/output) real ( kind = 8 ) array, dimension (LDAB,N)
!          On entry, the matrix A in band storage, in rows KL+1 to
!          2*KL+KU+1; rows 1 to KL of the array need not be set.
!          The j-th column of A is stored in the j-th column of the
!          array AB as follows:
!          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
!
!          On exit, details of the factorization: U is stored as an
!          upper triangular band matrix with KL+KU superdiagonals in
!          rows 1 to KL+KU+1, and the multipliers used during the
!          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
!          See below for further details.
!
!  LDAB    (input) integer
!          The leading dimension of the array AB.  2*KL+KU+1 <= LDAB.
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!
!  The band storage scheme is illustrated by the following example, when
!  M = N = 6, KL = 2, KU = 1:
!
!  On entry:                       On exit:
!
!      *    *    *    +    +    +       *    *    *   u14  u25  u36
!      *    *    +    +    +    +       *    *   u13  u24  u35  u46
!      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
!     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
!     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
!     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *
!
!  Array elements marked * are not used by the routine; elements marked
!  + need not be set on entry, but are required by the routine to store
!  elements of U because of fill-in resulting from the row interchanges.
!
  real ( kind = 8 )   one, zero
  parameter          ( one = 1.0d+0, zero = 0.0d+0 )
  integer ( kind = 4 )            nbmax, ldwork
  parameter          ( nbmax = 64, ldwork = nbmax+1 )
!
  integer ( kind = 4 )            i, i2, i3, ii, ip, j, j2, j3, jb, jj, jm, jp
  integer ( kind = 4 ) ju, k2, km, kv, nb, nw
  real ( kind = 8 )   temp
  real ( kind = 8 )   work13( ldwork, nbmax )
  real ( kind = 8 )   work31( ldwork, nbmax )
!
  integer ( kind = 4 )            idamax, ilaenv
  external           idamax, ilaenv
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
  kv = ku + kl
!
!     Test the input parameters.
!
  info = 0
  if ( m < 0 ) then
     info = -1
  else if ( n < 0 ) then
     info = -2
  else if ( kl < 0 ) then
     info = -3
  else if ( ku < 0 ) then
     info = -4
  else if ( ldab < kl+kv+1 ) then
     info = -6
  end if
  if ( info /= 0 ) then
     call xerbla( 'dgbtrf', -info )
     return
  end if
!
!  Quick return if possible
!
  if ( m == 0 .or. n.eq.0D+00 ) then
    return
  end if
!
!  Determine the block size for this environment
!
  nb = ilaenv( 1, 'dgbtrf', ' ', m, n, kl, ku )
!
!  The block size must not exceed the limit set by the size of the
!  local arrays WORK13 and WORK31.
!
  nb = min( nb, nbmax )

  if ( nb <= 1 .or. kl < nb ) then
!
!  Use unblocked code
!
     call dgbtf2( m, n, kl, ku, ab, ldab, ipiv, info )
  else
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
     do j = 1, nb
        do i = 1, j - 1
           work13( i, j ) = zero
        end do
     end do
!
!        Zero the subdiagonal elements of the work array WORK31
!
     do j = 1, nb
        do i = j + 1, nb
           work31( i, j ) = zero
        end do
     end do
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
     do j = ku + 2, min( kv, n )
        do i = kv - j + 2, kl
           ab( i, j ) = zero
        end do
     end do
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
     ju = 1
!
     do 180 j = 1, min( m, n ), nb
        jb = min( nb, min( m, n )-j+1 )
!
!           The active part of the matrix is partitioned
!
!              A11   A12   A13
!              A21   A22   A23
!              A31   A32   A33
!
!           Here A11, A21 and A31 denote the current block of JB columns
!           which is about to be factorized. The number of rows in the
!           partitioning are JB, I2, I3 respectively, and the numbers
!           of columns are JB, J2, J3. The superdiagonal elements of A13
!           and the subdiagonal elements of A31 lie outside the band.
!
        i2 = min( kl-jb, m-j-jb+1 )
        i3 = min( jb, m-j-kl+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
        do 80 jj = j, j + jb - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
           if ( jj+kv <= n ) then
              do i = 1, kl
                 ab( i, jj+kv ) = zero
              end do
           end if
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
           km = min( kl, m-jj )
           jp = idamax( km+1, ab( kv+1, jj ), 1 )
           ipiv( jj ) = jp + jj - j
           if ( ab( kv+jp, jj ) /= zero ) then
              ju = max( ju, min( jj+ku+jp-1, n ) )
              if ( jp /= 1 ) then
!
!                    Apply interchange to columns J to J+JB-1
!
                 if ( jp+jj-1 < j+kl ) then

                    call dswap( jb, ab( kv+1+jj-j, j ), ldab-1, &
                                   ab( kv+jp+jj-j, j ), ldab-1 )
                 else
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                    call dswap( jj-j, ab( kv+1+jj-j, j ), ldab-1, &
                                   work31( jp+jj-j-kl, 1 ), ldwork )
                    call dswap( j+jb-jj, ab( kv+1, jj ), ldab-1, &
                                   ab( kv+jp, jj ), ldab-1 )
                 end if
              end if
!
!  Compute multipliers
!
              call dscal( km, one / ab( kv+1, jj ), ab( kv+2, jj ), 1 )
!
!  Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
              jm = min( ju, j+jb-1 )

              if ( jj < jm ) then
                 call dger( km, jm-jj, -one, ab( kv+2, jj ), 1, &
                                ab( kv, jj+1 ), ldab-1, &
                                ab( kv+1, jj+1 ), ldab-1 )
              end if

           else
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
              if ( info == 0 ) then
                 info = jj
              end if

           end if
!
!  Copy current column of A31 into the work array WORK31
!
           nw = min( jj-j+1, i3 )
           if ( 0 < nw ) then
             call dcopy( nw, ab( kv+kl+1-jj+j, jj ), 1, &
                            work31( 1, jj-j+1 ), 1 )
           end if
   80       continue

        if ( j+jb <= n ) then
!
!              Apply the row interchanges to the other blocks.
!
           j2 = min( ju-j+1, kv ) - jb
           j3 = max( 0, ju-j-kv+1 )
!
!  Use DLASWP to apply the row interchanges to A12, A22, and A32.
!
           call dlaswp( j2, ab( kv+1-jb, j+jb ), ldab-1, 1, jb, ipiv( j ), 1 )
!
!              Adjust the pivot indices.
!
           do 90 i = j, j + jb - 1
              ipiv( i ) = ipiv( i ) + j - 1
   90          continue
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
           k2 = j - 1 + jb + j2
           do 110 i = 1, j3
              jj = k2 + i
              do 100 ii = j + i - 1, j + jb - 1
                 ip = ipiv( ii )
                 if ( ip /= ii ) then
                    temp = ab( kv+1+ii-jj, jj )
                    ab( kv+1+ii-jj, jj ) = ab( kv+1+ip-jj, jj )
                    ab( kv+1+ip-jj, jj ) = temp
                 end if
  100             continue
  110          continue
!
!  Update the relevant part of the trailing submatrix
!
           if ( 0 < j2 ) then
!
!                 Update A12
!
              call dtrsm( 'left', 'lower', 'no transpose', 'unit', &
                             jb, j2, one, ab( kv+1, j ), ldab-1, &
                             ab( kv+1-jb, j+jb ), ldab-1 )
!
              if ( 0 < i2 ) then
!
!                    Update A22
!
                 call dgemm( 'no transpose', 'no transpose', i2, j2, &
                                jb, -one, ab( kv+1+jb, j ), ldab-1, &
                                ab( kv+1-jb, j+jb ), ldab-1, one, &
                                ab( kv+1, j+jb ), ldab-1 )
              end if
!
              if ( 0 < i3 ) then
!
!                    Update A32
!
                 call dgemm( 'no transpose', 'no transpose', i3, j2, &
                                jb, -one, work31, ldwork, &
                                ab( kv+1-jb, j+jb ), ldab-1, one, &
                                ab( kv+kl+1-jb, j+jb ), ldab-1 )
              end if
           end if
!
           if ( 0 < j3 ) then
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
              do jj = 1, j3
                 do ii = jj, jb
                    work13( ii, jj ) = ab( ii-jj+1, jj+j+kv-1 )
                 end do
              end do
!
!                 Update A13 in the work array
!
              call dtrsm( 'left', 'lower', 'no transpose', 'unit', &
                             jb, j3, one, ab( kv+1, j ), ldab-1, &
                             work13, ldwork )
!
              if ( 0 < i2 ) then
!
!                    Update A23
!
                 call dgemm( 'no transpose', 'no transpose', i2, j3, &
                                jb, -one, ab( kv+1+jb, j ), ldab-1, &
                                work13, ldwork, one, ab( 1+jb, j+kv ), ldab-1 )
              end if

              if ( 0 < i3 ) then
!
!                    Update A33
!
                 call dgemm( 'no transpose', 'no transpose', i3, j3, &
                                jb, -one, work31, ldwork, work13, &
                                ldwork, one, ab( 1+kl, j+kv ), ldab-1 )
              end if
!
!                 Copy the lower triangle of A13 back into place
!
              do jj = 1, j3
                 do ii = jj, jb
                    ab( ii-jj+1, jj+j+kv-1 ) = work13( ii, jj )
                 end do
              end do
           end if
        else
!
!  Adjust the pivot indices.
!
           do i = j, j + jb - 1
              ipiv( i ) = ipiv( i ) + j - 1
           end do
        end if
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
        do 170 jj = j + jb - 1, j, -1
           jp = ipiv( jj ) - jj + 1
           if ( jp /= 1 ) then
!
!                 Apply interchange to columns J to JJ-1
!
              if ( jp+jj-1 < j+kl ) then
!
!                    The interchange does not affect A31
!
                 call dswap( jj-j, ab( kv+1+jj-j, j ), ldab-1, &
                               ab( kv+jp+jj-j, j ), ldab-1 )
              else
!
!  The interchange does affect A31
!
                 call dswap( jj-j, ab( kv+1+jj-j, j ), ldab-1, &
                             work31( jp+jj-j-kl, 1 ), ldwork )
              end if
           end if
!
!  Copy the current column of A31 back into place
!
           nw = min( i3, jj-j+1 )

           if ( 0 < nw ) then
              call dcopy( nw, work31( 1, jj-j+1 ), 1, &
                ab( kv+kl+1-jj+j, jj ), 1 )
           end if

  170       continue
  180    continue
  end if

  return
end
subroutine dgbtrs(trans,n,kl,ku,nrhs,ab,ldab,ipiv,b,ldb,info)
!
!***********************************************************************
!
!! DGBTRS ???
!
  character          trans
  integer ( kind = 4 )            info, kl, ku, ldab, ldb, n, nrhs
!
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   ab( ldab, * ), b( ldb, * )
!     ..
!  DGBTRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general band matrix A using the LU factorization computed
!  by DGBTRF.
!
!  Arguments
!
!  TRANS   (input) character*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) integer
!          The order of the matrix A.
!
!  KL      (input) integer
!          The number of subdiagonals within the band of A.
!
!  KU      (input) integer
!          The number of superdiagonals within the band of A.
!
!  NRHS    (input) integer
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.
!
!  AB      (input) real ( kind = 8 ) array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by DGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) integer
!          The leading dimension of the array AB.  2*KL+KU+1 <= LDAB.
!
!  IPIV    (input) integer array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) real ( kind = 8 ) array, dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.
!
!  LDB     (input) integer
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!
  real ( kind = 8 )   one
  parameter          ( one = 1.0d+0 )
!
  logical            lnoti, notran
  integer ( kind = 4 )            i, j, kd, l, lm
  logical            lsame
!
!     Test the input parameters.
!
  info = 0
  notran = lsame( trans, 'n' )
  if ( .not.notran .and. .not.lsame( trans, 't' ) .and. &
     .not. lsame( trans, 'c' ) ) then
     info = -1
  else if ( n < 0 ) then
     info = -2
  else if ( kl < 0 ) then
     info = -3
  else if ( ku < 0 ) then
     info = -4
  else if ( nrhs < 0 ) then
     info = -5
  else if ( ldab < ( 2*kl+ku+1 ) ) then
     info = -7
  else if ( ldb < max( 1, n ) ) then
     info = -10
  end if
  if ( info /= 0 ) then
     call xerbla( 'dgbtrs', -info )
     return
  end if
!
!  Quick return if possible
!
  if ( n == 0 .or. nrhs == 0 ) then
    return
  end if

  kd = ku + kl + 1
  lnoti = ( 0 < kl )

  if ( notran ) then
!
!  Solve  A*X = B.
!
!  Solve L*X = B, overwriting B with X.
!
!  L is represented as a product of permutations and unit lower
!  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!  where each transformation L(i) is a rank-one modification of
!  the identity matrix.
!
     if ( lnoti ) then
        do j = 1, n - 1
           lm = min( kl, n-j )
           l = ipiv( j )

           if ( l /= j ) then
              call dswap( nrhs, b( l, 1 ), ldb, b( j, 1 ), ldb )
           end if

           call dger( lm, nrhs, -one, ab( kd+1, j ), 1, b( j, 1 ), &
                         ldb, b( j+1, 1 ), ldb )
        end do
     end if

     do i = 1, nrhs
!
!  Solve U*X = B, overwriting B with X.
!
        call dtbsv( 'upper', 'no transpose', 'non-unit', n, kl+ku, &
          ab, ldab, b( 1, i ), 1 )
     end do

  else
!
!  Solve A'*X = B.
!
     do i = 1, nrhs
!
!  Solve U'*X = B, overwriting B with X.
!
        call dtbsv( 'upper', 'transpose', 'non-unit', n, kl+ku, ab, &
                       ldab, b( 1, i ), 1 )
     end do
!
!   Solve L'*X = B, overwriting B with X.
!
     if ( lnoti ) then
        do j = n - 1, 1, -1
           lm = min( kl, n-j )
           call dgemv( 'transpose', lm, nrhs, -one, b( j+1, 1 ), &
                          ldb, ab( kd+1, j ), 1, one, b( j, 1 ), ldb )
           l = ipiv( j )
           if ( l /= j ) then
             call dswap( nrhs, b( l, 1 ), ldb, b( j, 1 ), ldb )
           end if
        end do
     end if
  end if

  return
end
subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
!
!***********************************************************************
!
!! DGEMM ???
!
  character       TRANSA, TRANSB
  integer ( kind = 4 )            M, N, K, LDA, LDB, LDC
  real ( kind = 8 )   ALPHA, BETA
  real ( kind = 8 )   A( LDA, * ), B( LDB, * ), C( LDC, * )
!
!  DGEMM  performs one of the matrix-matrix operations
!
!     C : =  alpha*op( A )*op( B ) + beta*C,
!
!  where  op( X ) is one of
!
!     op( X ) = X   or   op( X ) = X',
!
!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!
!  Parameters
!
!  TRANSA - character*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n',  op( A ) = A.
!
!              TRANSA = 'T' or 't',  op( A ) = A'.
!
!              TRANSA = 'C' or 'c',  op( A ) = A'.
!
!           Unchanged on exit.
!
!  TRANSB - character*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSB = 'N' or 'n',  op( B ) = B.
!
!              TRANSB = 'T' or 't',  op( B ) = B'.
!
!              TRANSB = 'C' or 'c',  op( B ) = B'.
!
!           Unchanged on exit.
!
!  M      - integer.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.
!
!  K      - integer.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.
!
!  B      - real ( kind = 8 ) array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.
!
!  LDB    - integer.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.
!
!  BETA   - real ( kind = 8 ).
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.
!
!  C      - real ( kind = 8 ) array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).
!
!  LDC    - integer.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
!
  logical            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
  logical            NOTA, NOTB
  integer ( kind = 4 )            I, INFO, J, L, NROWA, NROWB
  real ( kind = 8 )   TEMP
!     .. Parameters ..
  real ( kind = 8 )   ONE         , zero
  PARAMETER        ( ONE = 1.0D+0, zero = 0.0D+0 )
!
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.
!
  NOTA  = LSAME( TRANSA, 'N' )
  NOTB  = LSAME( TRANSB, 'N' )
  if ( NOTA ) then
     NROWA = M
  else
     NROWA = K
  end if
  if ( NOTB ) then
     NROWB = K
  else
     NROWB = N
  end if
!
!     Test the input parameters.
!
  INFO = 0
  if ( ( .NOT. NOTA ) .and. &
       ( .NOT.LSAME( TRANSA, 'C' ) ) .and. &
       ( .NOT.LSAME( TRANSA, 'T' ) )      ) then
     INFO = 1
  else if ( ( .NOT.NOTB                 ).and. &
              ( .NOT.LSAME( TRANSB, 'C' ) ).and. &
              ( .NOT.LSAME( TRANSB, 'T' ) )      ) then
     INFO = 2
  else if ( M   < 0               ) then
     INFO = 3
  else if ( N   < 0               ) then
     INFO = 4
  else if ( K   < 0               ) then
     INFO = 5
  else if ( LDA < MAX( 1, NROWA ) ) then
     INFO = 8
  else if ( LDB < MAX( 1, NROWB ) ) then
     INFO = 10
  else if ( LDC < MAX( 1, M     ) ) then
     INFO = 13
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DGEMM ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if ( ( M == 0 ).OR.( N == 0 ).OR. &
     ( ( ( ALPHA == zero ).OR.( K == 0 ) ).and.( BETA == ONE ) ) ) then
    return
  end if
!
!     And if  alpha == zero.
!
  if ( ALPHA == zero ) then
     if ( BETA == zero ) then
        DO 20, J = 1, N
           DO 10, I = 1, M
              C( I, J ) = zero
   10          continue
   20       continue
     else
        DO 40, J = 1, N
           DO 30, I = 1, M
              C( I, J ) = BETA*C( I, J )
   30          continue
   40       continue
     end if
     return
  end if
!
!     Start the operations.
!
  if ( NOTB ) then
     if ( NOTA ) then
!
!           Form  C : =  alpha*A*B + beta*C.
!
        DO 90, J = 1, N
           if ( BETA == zero ) then
              DO 50, I = 1, M
                 C( I, J ) = zero
   50             continue
           else if ( BETA.NE.ONE ) then
              DO 60, I = 1, M
                 C( I, J ) = BETA*C( I, J )
   60             continue
           end if
           DO 80, L = 1, K
              if ( B( L, J ).NE.zero ) then
                 TEMP = ALPHA*B( L, J )
                 DO 70, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                continue
              end if
   80          continue
   90       continue
     else
!
!           Form  C : =  alpha*A'*B + beta*C
!
        DO 120, J = 1, N
           DO 110, I = 1, M
              TEMP = zero
              DO 100, L = 1, K
                 TEMP = TEMP + A( L, I )*B( L, J )
  100             continue
              if ( BETA == zero ) then
                 C( I, J ) = ALPHA*TEMP
              else
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
  110          continue
  120       continue
     end if
  else
     if ( NOTA ) then
!
!           Form  C : =  alpha*A*B' + beta*C
!
        DO 170, J = 1, N
           if ( BETA == zero ) then
              DO 130, I = 1, M
                 C( I, J ) = zero
  130             continue
           else if ( BETA.NE.ONE ) then
              DO 140, I = 1, M
                 C( I, J ) = BETA*C( I, J )
  140             continue
           end if
           DO 160, L = 1, K
              if ( B( J, L ).NE.zero ) then
                 TEMP = ALPHA*B( J, L )
                 DO 150, I = 1, M
                    C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                continue
              end if
  160          continue
  170       continue
     else
!
!           Form  C : =  alpha*A'*B' + beta*C
!
        DO J = 1, N
           DO I = 1, M
              TEMP = zero
              DO L = 1, K
                 TEMP = TEMP + A( L, I )*B( J, L )
              end do
              if ( BETA == zero ) then
                 C( I, J ) = ALPHA*TEMP
              else
                 C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
              end if
           end do
        end do
     end if
  end if

  return
end
subroutine dgemv( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
!
!***********************************************************************
!
!! DGEMV ???
!
!
  real ( kind = 8 )   ALPHA, BETA
  integer ( kind = 4 )            INCX, INCY, LDA, M, N
  character       TRANS
  real ( kind = 8 )   A( LDA, * ), X( * ), Y( * )
!
!  DGEMV  performs one of the matrix-vector operations
!
!     y : =  alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
!
!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.
!
!  Parameters
!
!  TRANS  - character*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!
!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
!
!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
!
!           Unchanged on exit.
!
!  M      - integer.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.
!
!  INCX   - integer.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  BETA   - real ( kind = 8 ).
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.
!
!  INCY   - integer.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!     .. Parameters ..
  real ( kind = 8 )   ONE         , zero
  PARAMETER        ( ONE = 1.0D+0, zero = 0.0D+0 )
!     .. Local Scalars ..
  real ( kind = 8 )   TEMP
  integer ( kind = 4 )            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
!     .. External Functions ..
  logical            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!
!     Test the input parameters.
!
  INFO = 0
  IF     ( .NOT.LSAME( TRANS, 'N' ).and. &
              .NOT.LSAME( TRANS, 'T' ).and. &
              .NOT.LSAME( TRANS, 'C' )      ) then
     INFO = 1
  else if ( M < 0 ) then
     INFO = 2
  else if ( N < 0 ) then
     INFO = 3
  else if ( LDA < MAX( 1, M ) ) then
     INFO = 6
  else if ( INCX == 0 ) then
     INFO = 8
  else if ( INCY == 0 ) then
     INFO = 11
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DGEMV ', INFO )
     return
  end if

  if ( ( M == 0 ).OR.( N == 0 ).OR. &
       ( ( ALPHA == zero ).and.( BETA == ONE ) ) ) then
    return
  end if
!
!  Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!  up the start points in  X  and  Y.
!
  if ( LSAME( TRANS, 'N' ) ) then
     LENX = N
     LENY = M
  else
     LENX = M
     LENY = N
  end if
  if ( 0 < INCX ) then
     KX = 1
  else
     KX = 1 - ( LENX - 1 )*INCX
  end if
  if ( 0 < INCY ) then
     KY = 1
  else
     KY = 1 - ( LENY - 1 )*INCY
  end if
!
!  Start the operations. In this version the elements of A are
!  accessed sequentially with one pass through A.
!
!  First form  y : =  beta*y.
!
  if ( BETA.NE.ONE ) then
     if ( INCY == 1 ) then
        if ( BETA == zero ) then
           DO I = 1, LENY
              Y( I ) = zero
           end do
        else
           DO 20, I = 1, LENY
              Y( I ) = BETA*Y( I )
   20          continue
        end if
     else
        IY = KY
        if ( BETA == zero ) then
           DO 30, I = 1, LENY
              Y( IY ) = zero
              IY      = IY   + INCY
   30          continue
        else
           DO 40, I = 1, LENY
              Y( IY ) = BETA*Y( IY )
              IY      = IY           + INCY
   40          continue
        end if
     end if
  end if

  if ( ALPHA == zero ) then
    return
  end if

  if ( LSAME( TRANS, 'N' ) ) then
!
!  Form  y : =  alpha*A*x + y.
!
     JX = KX
     if ( INCY == 1 ) then
        DO 60, J = 1, N
           if ( X( JX ).NE.zero ) then
              TEMP = ALPHA*X( JX )
              DO 50, I = 1, M
                 Y( I ) = Y( I ) + TEMP*A( I, J )
   50             continue
           end if
           JX = JX + INCX
   60       continue
     else
        DO 80, J = 1, N
           if ( X( JX ).NE.zero ) then
              TEMP = ALPHA*X( JX )
              IY   = KY
              DO 70, I = 1, M
                 Y( IY ) = Y( IY ) + TEMP*A( I, J )
                 IY      = IY      + INCY
   70             continue
           end if
           JX = JX + INCX
   80       continue
     end if
  else
!
!  Form  y : =  alpha*A'*x + y.
!
     JY = KY
     if ( INCX == 1 ) then
        DO 100, J = 1, N
           TEMP = zero
           DO 90, I = 1, M
              TEMP = TEMP + A( I, J )*X( I )
   90          continue
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
  100       continue
     else
        DO 120, J = 1, N
           TEMP = zero
           IX   = KX
           DO I = 1, M
              TEMP = TEMP + A( I, J )*X( IX )
              IX   = IX   + INCX
           end do
           Y( JY ) = Y( JY ) + ALPHA*TEMP
           JY      = JY      + INCY
  120       continue
     end if
  end if

  return
end
subroutine dgeqr2(m,n,a,lda,tau,work,info)
!
!***********************************************************************
!
!! DGEQR2 computes a QR factorization of a real m by n matrix A:
!  A = Q * R.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix A. 
!
!  N       (input) integer
!          The number of columns of the matrix A. 
!
!  A       (input/output) real ( kind = 8 ) array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(m,n) by n upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) real ( kind = 8 ) array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) real ( kind = 8 ) array, dimension (N)
!
!  INFO    (output) integer
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!
  real ( kind = 8 ) one
  parameter        ( one = 1.0D+0 )
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) a(lda,n)
  real ( kind = 8 ) aii
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real ( kind = 8 ) tau(*)
  real ( kind = 8 ) work(n)
  EXTERNAL           DLARF, DLARFG, XERBLA
!
!  Test the input arguments
!
  info = 0

  if ( m < 0) then
     info = -1
  else if ( n < 0) then
     info = -2
  else if ( LDA < MAX( 1, M ) ) then
     INFO = -4
  end if

  if ( info /= 0 ) then
     CALL XERBLA( 'DGEQR2', -INFO )
     return
  end if

  k = MIN( M, N )
!
!  Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
  DO I = 1, K

     call dlarfg(m-i+1,a(i,i),a(min(i+1, m),i),1,tau(i))
!
!  Apply H(i) to A(i:m,i+1:n) from the left
!
     if ( i < n) then
        aii = a( i,i )
        a( i,i ) = one
        call dlarf( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
          A( I, I+1 ), LDA, WORK )
        A( I, I ) = AII
     end if

  end do

  return
end
subroutine dgeqrf( m, n, a, lda, tau, work, lwork, info )
!
!***********************************************************************
!
!! DGEQRF computes a QR factorization of a real M-by-N matrix A:
!  A = Q * R.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix A.
!
!  N       (input) integer
!          The number of columns of the matrix A.
!
!  A       (input/output) real ( kind = 8 ) array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of min(m,n) elementary reflectors (see Further
!          Details).
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) real ( kind = 8 ) array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) real ( kind = 8 ) array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) integer
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is
!          the optimal blocksize.
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!
!  The matrix Q is represented as a product of elementary reflectors
!
!     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!  and tau in TAU(i).
!
!
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) lwork
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
!
  real ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) iinfo
  integer ( kind = 4 ) ilaenv
  integer ( kind = 4 ) info
  integer ( kind = 4 ) iws
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ldwork
  integer ( kind = 4 ) nb
  integer ( kind = 4 ) nbmin
  integer ( kind = 4 ) nx
  real ( kind = 8 ) tau(min(m,n))
  real ( kind = 8 ) work(lwork)
!
  EXTERNAL           DGEQR2, DLARFB, DLARFT, XERBLA
  EXTERNAL           ILAENV
!
!  Test the input arguments
!
  info = 0

  if ( M < 0 ) then
     INFO = -1
  else if ( N < 0 ) then
     INFO = -2
  else if ( LDA < MAX( 1, M ) ) then
     INFO = -4
  else if ( LWORK < MAX( 1, N ) ) then
     INFO = -7
  end if

  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DGEQRF', -INFO )
     return
  end if
!
!  Quick return if possible
!
  K = MIN( M, N )

  if ( K == 0 ) then
     WORK( 1 ) = 1
     return
  end if
!
!  Determine the block size.
!
  NB = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
  nbmin = 2
  NX = 0
  IWS = N
!
!  Determine when to cross over from blocked to unblocked code.
!
  if ( NB > 1 .and. NB < K ) then
     NX = MAX( 0, ILAENV( 3, 'DGEQRF', ' ', M, N, -1, -1 ) )
!
!  Determine if workspace is large enough for blocked code.
!
     if ( nx < k) then
        LDWORK = N
        IWS = LDWORK*NB
!
!  Not enough workspace to use optimal NB:  reduce NB and
!  determine the minimum value of NB.
!
        if ( LWORK < IWS ) then
           NB = LWORK / LDWORK
           NBMIN = MAX(2,ILAENV( 2,'DGEQRF',' ', M, N, -1,-1))
        end if

     end if

  end if

!
!  Use blocked code initially
!
  if ( nb >= nbmin.and.nb < k.and.nx < k) then
     do i = 1, k - nx, nb
        ib = min( k-i+1, nb )
!
!  Compute the QR factorization of the current block A(i:m,i:i+ib-1)
!
        call dgeqr2(m-i+1,ib,a(i,i),lda,tau(i),work,iinfo)
!
!  Form the triangular factor of the block reflector
!  H = H(i) H(i+1) . . . H(i+ib-1)
!
        if ( i+ib <= n) then
           call dlarft( 'Forward', 'Columnwise', M-I+1, IB, &
                      a( i,i ), lda, tau( i ), work, ldwork )
!
!  Apply H' to A(i:m,i+ib:n) from the left
!
           call dlarfb( 'Left', 'Transpose', 'Forward', &
             'Columnwise', m-i+1, n-i-ib+1, ib, &
             a( i,i ), lda, work, ldwork, a( i,i+ib ), &
             lda, work( ib+1 ), ldwork )
        end if
    end do

  else

     i = 1

  end if

!
!  Use unblocked code to factor the last or only block.
!
  if ( i <= k) then
    call dgeqr2(m-i+1,n-i+1,a(i,i),lda,tau(i),work,iinfo)
  end if

  work(1) = iws

  return
end
subroutine dger( m, n, alpha, x, incx, y, incy, a, lda )
!
!***********************************************************************
!
!! DGER ???
!
  real ( kind = 8 )   ALPHA
  integer ( kind = 4 )            INCX, INCY, LDA, M, N
  real ( kind = 8 )   A( LDA, * ), X( * ), Y( * )
!
!  DGER   performs the rank 1 operation
!
!     A : =  alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!
!  M      - integer.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - integer.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - integer.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!     .. Parameters ..
  real ( kind = 8 )   zero
  PARAMETER        ( zero = 0.0D+0 )
!     .. Local Scalars ..
  real ( kind = 8 )   TEMP
  integer ( kind = 4 )            I, INFO, IX, J, JY, KX
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!
!
!     Test the input parameters.
!
  INFO = 0
  IF     ( M < 0 ) then
     INFO = 1
  else if ( N < 0 ) then
     INFO = 2
     write ( *, * ) ' '
     write ( *, * ) 'DGER - Fatal error!'
     write ( *, * ) '  Input value of N was less than 0.'
     write ( *, * ) '  N = ',n
  else if ( incx == 0 ) then
     INFO = 5
  else if ( INCY == 0 ) then
     INFO = 7
  else if ( LDA < MAX( 1, M ) ) then
     INFO = 9
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DGER  ', INFO )
     return
  end if

  if ( ( M == 0 ).OR.( N == 0 ).OR.( ALPHA == zero ) ) then
    return
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  if ( INCY > 0 ) then
     JY = 1
  else
     JY = 1 - ( N - 1 )*INCY
  end if
  if ( INCX == 1 ) then
     DO 20, J = 1, N
        if ( Y( JY ).NE.zero ) then
           TEMP = ALPHA*Y( JY )
           DO 10, I = 1, M
              A( I, J ) = A( I, J ) + X( I )*TEMP
   10          continue
        end if
        JY = JY + INCY
   20    continue
  else
     if ( INCX > 0 ) then
        KX = 1
     else
        KX = 1 - ( M - 1 )*INCX
     end if
     DO 40, J = 1, N
        if ( Y( JY ).NE.zero ) then
           TEMP = ALPHA*Y( JY )
           IX   = KX
           DO I = 1, M
              A( I, J ) = A( I, J ) + X( IX )*TEMP
              IX        = IX        + INCX
           end do
        end if
        JY = JY + INCY
   40    continue
  end if

  return
end
subroutine dgetf2(m,n,a,lda,ipiv,info)
!
!***********************************************************************
!
!! DGETF2 ???
!
!
  integer ( kind = 4 )            INFO, LDA, M, N
  integer ( kind = 4 )            IPIV( * )
  real ( kind = 8 )   A( LDA, * )
!
!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 2 BLAS version of the algorithm.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real ( kind = 8 ) array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
  real ( kind = 8 )   ONE, zero
  PARAMETER          ( ONE = 1.0D+0, zero = 0.0D+0 )
!
  integer ( kind = 4 )            J, JP
!     ..
!     .. External Functions ..
  integer ( kind = 4 )            IDAMAX
  EXTERNAL           IDAMAX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
!
!     Test the input parameters.
!
  INFO = 0
  if ( M < 0 ) then
     INFO = -1
  else if ( N < 0 ) then
     INFO = -2
  else if ( LDA < MAX( 1, M ) ) then
     INFO = -4
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DGETF2', -INFO )
     return
  end if

  if ( M == 0 .OR. N == 0 ) then
    return
  end if
!
  DO 10 J = 1, MIN( M, N )
!
!        Find pivot and test for singularity.
!
     JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
     IPIV( J ) = JP
     if ( A( JP, J ).NE.zero ) then
!
!  Apply the interchange to columns 1:N.
!
        if ( JP.NE.J ) then
           CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
        end if
!
!   Compute elements J+1:M of J-th column.
!
        if ( J < M ) then
           CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
        end if
     else if ( INFO == 0 ) then
        INFO = J
     end if

     if ( J < MIN( M, N ) ) then
!
!  Update trailing submatrix.
!
        call dger( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA, &
          A( J+1, J+1 ), LDA )

     end if

   10 continue

  return
end
subroutine dgetrf( M, N, A, LDA, IPIV, INFO )
!
!***********************************************************************
!
!! DGETRF ???
!
!
  integer ( kind = 4 )            INFO, LDA, M, N
  integer ( kind = 4 )            IPIV( * )
  real ( kind = 8 )   A( LDA, * )
!
!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) real ( kind = 8 ) array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) integer array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
  real ( kind = 8 )   ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  integer ( kind = 4 )            I, IINFO, J, JB, NB
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
!     ..
!     .. External Functions ..
  integer ( kind = 4 )            ILAENV
  EXTERNAL           ILAENV
!
!     Test the input parameters.
!
  INFO = 0
  if ( M < 0 ) then
     INFO = -1
  else if ( N < 0 ) then
     INFO = -2
  else if ( LDA < MAX( 1, M ) ) then
     INFO = -4
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DGETRF', -INFO )
     return
  end if

  if ( M == 0 .OR. N == 0 ) then
    return
  end if
!
!     Determine the block size for this environment.
!
  NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
  if ( NB <= 1 .OR. NB >= MIN( M, N ) ) then
!
!        Use unblocked code.
!
     CALL DGETF2( M, N, A, LDA, IPIV, INFO )
  else
!
!        Use blocked code.
!
     DO 20 J = 1, MIN( M, N ), NB
        JB = MIN( MIN( M, N )-J+1, NB )
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
        CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
!
!           Adjust INFO and the pivot indices.
!
        if ( INFO == 0 .and. IINFO > 0 ) then
          INFO = IINFO + J - 1
        end if

        DO I = J, MIN( M, J+JB-1 )
           IPIV( I ) = J - 1 + IPIV( I )
        end do
!
!           Apply interchanges to columns 1:J-1.
!
        CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )

        if ( J+JB <= N ) then
!
!  Apply interchanges to columns J+JB:N.
!
           CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1, IPIV, 1 )
!
!  Compute block row of U.
!
           CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,&
             N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ), LDA )
           if ( J+JB <= M ) then
!
!  Update trailing submatrix.
!
              CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1, &
                N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA, &
                A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ), LDA )
           end if
        end if
   20    continue
  end if

  return
end
subroutine dgetrs( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!***********************************************************************
!
!! DGETRS ???
!
  character          TRANS
  integer ( kind = 4 )            INFO, LDA, LDB, N, NRHS
  integer ( kind = 4 )            IPIV( * )
  real ( kind = 8 )   A( LDA, * ), B( LDB, * )
!
!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general N-by-N matrix A using the LU factorization computed
!  by DGETRF.
!
!  Arguments
!
!  TRANS   (input) character*1
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) integer
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) integer
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  A       (input) real ( kind = 8 ) array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.
!
!  LDA     (input) integer
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) integer array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  B       (input/output) real ( kind = 8 ) array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) integer
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
  real ( kind = 8 )   ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  logical            NOTRAN
!     ..
!     .. External Functions ..
  logical            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLASWP, DTRSM, XERBLA
!
!     Test the input parameters.
!
  INFO = 0
  NOTRAN = LSAME( TRANS, 'N' )
  if ( .NOT.NOTRAN .and. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
    LSAME( TRANS, 'C' ) ) then
     INFO = -1
  else if ( N < 0 ) then
     INFO = -2
  else if ( NRHS < 0 ) then
     INFO = -3
  else if ( LDA < MAX( 1, N ) ) then
     INFO = -5
  else if ( LDB < MAX( 1, N ) ) then
     INFO = -8
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DGETRS', -INFO )
     return
  end if
!
!     Quick return if possible
!
  if ( N == 0 .OR. NRHS == 0 ) then
    return
  end if

  if ( NOTRAN ) then
!
!        Solve A * X = B.
!
!        Apply row interchanges to the right hand sides.
!
     CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
!
!        Solve L*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS, &
       ONE, A, LDA, B, LDB )
!
!        Solve U*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N, &
                  NRHS, ONE, A, LDA, B, LDB )
  else
!
!        Solve A' * X = B.
!
!        Solve U'*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS, &
       ONE, A, LDA, B, LDB )
!
!        Solve L'*X = B, overwriting B with X.
!
     CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE, &
       A, LDA, B, LDB )
!
!        Apply row interchanges to the solution vectors.
!
     CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
  end if

  return
end
function dlapy2( x, y )
!
!***********************************************************************
!
!! DLAPY2 returns sqrt(x**2+y**2), avoiding overflow.
!
!
!  Arguments
!
!  X       (input) real ( kind = 8 )
!  Y       (input) real ( kind = 8 )
!          X and Y specify the values x and y.
!
  real ( kind = 8 ) dlapy2
  real ( kind = 8 )   X, Y
  real ( kind = 8 )   zero
  PARAMETER          ( zero = 0.0D0 )
  real ( kind = 8 )   ONE
  PARAMETER          ( ONE = 1.0D0 )
  real ( kind = 8 )   W, XABS, YABS, Z
!
  XABS = ABS( X )
  YABS = ABS( Y )
  W = MAX( XABS, YABS )
  Z = MIN( XABS, YABS )

  if ( Z == 0.0D+00 ) then
    DLAPY2 = W
  else
    DLAPY2 = W*SQRT( ONE+( Z / W )**2 )
  end if

  return
end
subroutine dlarf( side, m, n, v, incv, tau, c, ldc, work )
!
!***********************************************************************
!
!! DLARF ???
!
  character          SIDE
  integer ( kind = 4 )            INCV, LDC, M, N
  real ( kind = 8 )   TAU
  real ( kind = 8 )   C( LDC, * ), V( * ), WORK( * )
!
!
!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!
!  SIDE    (input) character*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) integer
!          The number of rows of the matrix C.
!
!  N       (input) integer
!          The number of columns of the matrix C.
!
!  V       (input) real ( kind = 8 ) array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) integer
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) real ( kind = 8 )
!          The value tau in the representation of H.
!
!  C       (input/output) real ( kind = 8 ) array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) integer
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) real ( kind = 8 ) array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!
!     .. Parameters ..
  real ( kind = 8 )   ONE, zero
  PARAMETER          ( ONE = 1.0D+0, zero = 0.0D+0 )
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMV, DGER
!     ..
!     .. External Functions ..
  logical            LSAME
  EXTERNAL           LSAME
!
  if ( LSAME( SIDE, 'L' ) ) then
!
!  Form  H * C
!
     if ( TAU /= 0.0D+00 ) then
!
!  w : =  C' * v
!
        CALL DGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, zero, WORK, 1 )
!
!  C : =  C - v * w'
!
        CALL DGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
     end if
  else
!
!  Form  C * H
!
     if ( TAU /= 0.0D+00 ) then
!
!           w : =  C * v
!
        CALL DGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, zero, WORK, 1 )
!
!  C : =  C - w * v'
!
        CALL DGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
     end if
  end if

  return
end
subroutine dlarfb( side, trans, direct, storev, m, n, k, v, ldv, &
  t, ldt, c, ldc, work, ldwork )
!
!***********************************************************************
!
!! DLARFB ???
!
  character          DIRECT, SIDE, STOREV, TRANS
  integer ( kind = 4 )            K, LDC, LDT, LDV, LDWORK, M, N
  real ( kind = 8 )   C( LDC, * ), T( LDT, * ), V( LDV, * )
  real ( kind = 8 ) WORK( LDWORK, * )
!
!  DLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.
!
!  Arguments
!
!  SIDE    (input) character*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right
!
!  TRANS   (input) character*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)
!
!  DIRECT  (input) character*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) character*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise
!
!  M       (input) integer
!          The number of rows of the matrix C.
!
!  N       (input) integer
!          The number of columns of the matrix C.
!
!  K       (input) integer
!          The order of the matrix T ( =  the number of elementary
!          reflectors whose product defines the block reflector).
!
!  V       (input) real ( kind = 8 ) array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) integer
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.
!
!  T       (input) real ( kind = 8 ) array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.
!
!  LDT     (input) integer
!          The leading dimension of the array T. LDT >= K.
!
!  C       (input/output) real ( kind = 8 ) array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
!
!  LDC     (input) integer
!          The leading dimension of the array C. LDA >= max(1,M).
!
!  WORK    (workspace) real ( kind = 8 ) array, dimension (LDWORK,K)
!
!  LDWORK  (input) integer
!          The leading dimension of the array WORK.
!          If SIDE = 'L', LDWORK >= max(1,N);
!          if SIDE = 'R', LDWORK >= max(1,M).
!
  real ( kind = 8 )   ONE
  PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
  character          TRANST
  integer ( kind = 4 )            I, J
!
  logical            LSAME
  EXTERNAL           LSAME
!
  EXTERNAL           DCOPY, DGEMM, DTRMM
!
!     Quick return if possible
!
  if ( M <= 0 .OR. N <= 0 ) then
    return
  end if

  if ( LSAME( TRANS, 'N' ) ) then
     TRANST = 'T'
  else
     TRANST = 'N'
  end if

  if ( LSAME( STOREV, 'C' ) ) then
!
     if ( LSAME( DIRECT, 'F' ) ) then
!
!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.
!
        if ( LSAME( SIDE, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W : =  C' * V = (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W : =  C1'
!
           DO J = 1, K
             CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
           end do
!
!  W : =  W * V1
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
             K, ONE, V, LDV, WORK, LDWORK )

           if ( M > K ) then
!
!                 W : =  W + C2'*V2
!
              CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                            ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
             ONE, T, LDT, WORK, LDWORK )
!
!  C : =  C - V * W'
!
           if ( M > K ) then
!
!  C2 : =  C2 - V2 * W'
!
              CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, C( K+1, 1 ), LDC )
           end if
!
!              W : =  W * V1'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
             ONE, V, LDV, WORK, LDWORK )
!
!              C1 : =  C1 - W'
!
           DO J = 1, K
              DO I = 1, N
                 C( J, I ) = C( J, I ) - WORK( I, J )
              end do
           end do

        else if ( LSAME( SIDE, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W : =  C * V = (C1*V1 + C2*V2)  (stored in WORK)
!
!              W : =  C1
!
           DO J = 1, K
              CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
           end do
!
!  W : =  W * V1
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
             K, ONE, V, LDV, WORK, LDWORK )

           if ( N > K ) then
!
!                 W : =  W + C2 * V2
!
              CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                            ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, &
                            ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
             ONE, T, LDT, WORK, LDWORK )
!
!  C : =  C - W * V'
!
           if ( N > K ) then
!
!                 C2 : =  C2 - W * V2'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                             -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, &
                             C( 1, K+1 ), LDC )
           end if
!
!              W : =  W * V1'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                          ONE, V, LDV, WORK, LDWORK )
!
!              C1 : =  C1 - W
!
           DO J = 1, K
              DO I = 1, M
                 C( I, J ) = C( I, J ) - WORK( I, J )
              end do
           end do
        end if

     else
!
!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.
!
        if ( LSAME( SIDE, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W : =  C' * V = (C1'*V1 + C2'*V2)  (stored in WORK)
!
!              W : =  C2'
!
           DO J = 1, K
              CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
           end do
!
!              W : =  W * V2
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
             K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )

           if ( M > K ) then
!
!                 W : =  W + C1'*V1
!
              CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K, &
                             ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                         ONE, T, LDT, WORK, LDWORK )
!
!              C : =  C - V * W'
!
           if ( M > K ) then
!
!                 C1 : =  C1 - V1 * W'
!
              CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K, &
                -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC )
           end if
!
!              W : =  W * V2'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                          ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 : =  C2 - W'
!
           DO 90 J = 1, K
              DO 80 I = 1, N
                 C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
   80             continue
   90          continue
!
        else if ( LSAME( SIDE, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W : =  C * V = (C1*V1 + C2*V2)  (stored in WORK)
!
!              W : =  C2
!
           DO J = 1, K
              CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
           end do
!
!              W : =  W * V2
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
             K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
           if ( N > K ) then
!
!                 W : =  W + C1 * V1
!
              CALL DGEMM( 'No transpose', 'No transpose', M, K, N-K, &
                ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
             ONE, T, LDT, WORK, LDWORK )
!
!              C : =  C - W * V'
!
           if ( N > K ) then
!
!                 C1 : =  C1 - W * V1'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, N-K, K, &
                -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
           end if
!
!              W : =  W * V2'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
             ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
!
!              C2 : =  C2 - W
!
           DO 120 J = 1, K
              DO 110 I = 1, M
                 C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             continue
  120          continue
        end if
     end if
!
  else if ( LSAME( STOREV, 'R' ) ) then
!
     if ( LSAME( DIRECT, 'F' ) ) then
!
!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.
!
        if ( LSAME( SIDE, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W : =  C' * V' = (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W : =  C1'
!
           DO 130 J = 1, K
              CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
  130          continue
!
!              W : =  W * V1'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', N, K, &
                    ONE, V, LDV, WORK, LDWORK )
           if ( M > K ) then
!
!                 W : =  W + C2'*V2'
!
              CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K, &
                 ONE, T, LDT, WORK, LDWORK )
!
!              C : =  C - V' * W'
!
           if ( M > K ) then
!
!                 C2 : =  C2 - V2' * W'
!
              CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                V( 1, K+1 ), LDV, WORK, LDWORK, ONE, C( K+1, 1 ), LDC )
           end if
!
!              W : =  W * V1
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N, &
                  K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 : =  C1 - W'
!
           DO 150 J = 1, K
              DO 140 I = 1, N
                 C( J, I ) = C( J, I ) - WORK( I, J )
  140             continue
  150          continue
!
        else if ( LSAME( SIDE, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W : =  C * V' = (C1*V1' + C2*V2')  (stored in WORK)
!
!              W : =  C1
!
           DO J = 1, K
              CALL DCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
           end do
!
!              W : =  W * V1'
!
           CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Unit', M, K, &
                      ONE, V, LDV, WORK, LDWORK )
           if ( N > K ) then
!
!                 W : =  W + C2 * V2'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K, &
             ONE, T, LDT, WORK, LDWORK )
!
!              C : =  C - W * V
!
           if ( N > K ) then
!
!                 C2 : =  C2 - W * V2
!
              CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, C( 1, K+1 ), LDC )
           end if
!
!              W : =  W * V1
!
           CALL DTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M, &
                         K, ONE, V, LDV, WORK, LDWORK )
!
!              C1 : =  C1 - W
!
           DO J = 1, K
              DO I = 1, M
                 C( I, J ) = C( I, J ) - WORK( I, J )
              end do
           end do

        end if

     else
!
!  Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.
!
        if ( LSAME( SIDE, 'L' ) ) then
!
!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )
!
!              W : =  C' * V' = (C1'*V1' + C2'*V2') (stored in WORK)
!
!              W : =  C2'
!
           DO 190 J = 1, K
              CALL DCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
  190          continue
!
!              W : =  W * V2'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K, &
             ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )

           if ( M > K ) then
!
!                 W : =  W + C1'*V1'
!
              CALL DGEMM( 'Transpose', 'Transpose', N, K, M-K, ONE, &
                         C, LDC, V, LDV, ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T'  or  W * T
!
           CALL DTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K, &
                       ONE, T, LDT, WORK, LDWORK )
!
!              C : =  C - V' * W'
!
           if ( M > K ) then
!
!                 C1 : =  C1 - V1' * W'
!
              CALL DGEMM( 'Transpose', 'Transpose', M-K, N, K, -ONE, &
                V, LDV, WORK, LDWORK, ONE, C, LDC )
           end if
!
!              W : =  W * V2
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N, &
                       K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
!
!              C2 : =  C2 - W'
!
           DO J = 1, K
              DO I = 1, N
                 C( M-K+J, I ) = C( M-K+J, I ) - WORK( I, J )
              end do
           end do

        else if ( LSAME( SIDE, 'R' ) ) then
!
!              Form  C * H  or  C * H'  where  C = ( C1  C2 )
!
!              W : =  C * V' = (C1*V1' + C2*V2')  (stored in WORK)
!
!              W : =  C2
!
           DO J = 1, K
              CALL DCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
           end do
!
!              W : =  W * V2'
!
           CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', M, K, &
                          ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
           if ( N > K ) then
!
!                 W : =  W + C1 * V1'
!
              CALL DGEMM( 'No transpose', 'Transpose', M, K, N-K, &
                             ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
           end if
!
!              W : =  W * T  or  W * T'
!
           CALL DTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K, &
             ONE, T, LDT, WORK, LDWORK )
!
!              C : =  C - W * V
!
           if ( N > K ) then
!
!                 C1 : =  C1 - W * V1
!
              CALL DGEMM( 'No transpose', 'No transpose', M, N-K, K, &
                            -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
           end if
!
!              W : =  W * V2
!
           CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M, &
             K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
!
!              C1 : =  C1 - W
!
           DO J = 1, K
              DO I = 1, M
                 C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
              end do
           end do

        end if

     end if
  end if

  return
end
subroutine dlarft( direct, storev, n, k, v, ldv, tau, t, ldt )
!
!***********************************************************************
!
!! DLARFT ???
!
  character          DIRECT, STOREV
  integer ( kind = 4 )            K, LDT, LDV, N
  real ( kind = 8 )   T( LDT, * ), TAU( * ), V( LDV, * )
!
!  DLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.
!
!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
!
!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
!
!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and
!
!     H  =  I - V * T * V'
!
!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and
!
!     H  =  I - V' * T * V
!
!  Arguments
!
!  DIRECT  (input) character*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)
!
!  STOREV  (input) character*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise
!
!  N       (input) integer
!          The order of the block reflector H. N >= 0.
!
!  K       (input) integer
!          The order of the triangular factor T ( =  the number of
!          elementary reflectors). K >= 1.
!
!  V       (input/output) real ( kind = 8 ) array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.
!
!  LDV     (input) integer
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
!
!  TAU     (input) real ( kind = 8 ) array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).
!
!  T       (output) real ( kind = 8 ) array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.
!
!  LDT     (input) integer
!          The leading dimension of the array T. LDT >= K.
!
!  Further Details
!
!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.
!
!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
!
!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )
!
!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
!
!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )
!
  real ( kind = 8 )   ONE, zero
  PARAMETER          ( ONE = 1.0D+0, zero = 0.0D+0 )
!
  integer ( kind = 4 )            I, J
  real ( kind = 8 )   VII
!     ..
!     .. External Subroutines ..
  EXTERNAL           DGEMV, DTRMV
!     ..
!     .. External Functions ..
  logical            LSAME
  EXTERNAL           LSAME
!
  if ( N == 0 ) then
    return
  end if

  if ( LSAME( DIRECT, 'F' ) ) then
     DO 20 I = 1, K
        if ( TAU( I ) == zero ) then
!
!              H(i)  =  I
!
           DO J = 1, I
              T( J, I ) = zero
           end do
        else
!
!              general case
!
           VII = V( I, I )
           V( I, I ) = ONE
           if ( LSAME( STOREV, 'C' ) ) then
!
!                 T(1:i-1,i) : =  - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
!
              CALL DGEMV( 'Transpose', N-I+1, I-1, -TAU( I ), &
                V( I, 1 ), LDV, V( I, I ), 1, zero, T( 1, I ), 1 )
           else
!
!                 T(1:i-1,i) : =  - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
!
              CALL DGEMV( 'No transpose', I-1, N-I+1, -TAU( I ), &
                V( 1, I ), LDV, V( I, I ), LDV, zero, T( 1, I ), 1 )
           end if
           V( I, I ) = VII
!
!              T(1:i-1,i) : =  T(1:i-1,1:i-1) * T(1:i-1,i)
!
           CALL DTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, &
             LDT, T( 1, I ), 1 )
           T( I, I ) = TAU( I )
        end if
   20    continue
  else
     DO 40 I = K, 1, -1
        if ( TAU( I ) == zero ) then
!
!              H(i)  =  I
!
           DO J = I, K
              T( J, I ) = zero
           end do
        else
!
!              general case
!
           if ( I < K ) then
              if ( LSAME( STOREV, 'C' ) ) then
                 VII = V( N-K+I, I )
                 V( N-K+I, I ) = ONE
!
!                    T(i+1:k,i) : =
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
!
                 CALL DGEMV( 'Transpose', N-K+I, K-I, -TAU( I ), &
                                V( 1, I+1 ), LDV, V( 1, I ), 1, zero, &
                               T( I+1, I ), 1 )
                 V( N-K+I, I ) = VII
              else
                 VII = V( I, N-K+I )
                 V( I, N-K+I ) = ONE
!
!                    T(i+1:k,i) : =
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
!
                 CALL DGEMV( 'No transpose', K-I, N-K+I, -TAU( I ), &
                   V( I+1, 1 ), LDV, V( I, 1 ), LDV, zero, T( I+1, I ), 1 )
                 V( I, N-K+I ) = VII
              end if
!
!                 T(i+1:k,i) : =  T(i+1:k,i+1:k) * T(i+1:k,i)
!
              CALL DTRMV( 'Lower', 'No transpose', 'Non-unit', K-I, &
                            T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
           end if
           T( I, I ) = TAU( I )
        end if
   40    continue
  end if

  return
end
subroutine dlaswp(n,a,lda,k1,k2,ipiv,incx)
!
!***********************************************************************
!
!! DLASWP ???
!
  integer ( kind = 4 )            incx, k1, k2, lda, n
  integer ( kind = 4 )            ipiv( * )
  real ( kind = 8 )   a( lda, * )
!
!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!
!  N       (input) integer
!          The number of columns of the matrix A.
!
!  A       (input/output) real ( kind = 8 ) array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) integer
!          The leading dimension of the array A.
!
!  K1      (input) integer
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) integer
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) integer array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) integer
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
!
!
  integer ( kind = 4 )            i, ip, ix
!
  external           dswap
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
  if ( incx == 0 ) then
    return
  end if
  if ( incx > 0 ) then
     ix = k1
  else
     ix = 1 + ( 1-k2 )*incx
  end if
  if ( incx == 1 ) then
     do i = k1, k2
        ip = ipiv( i )
        if ( ip /= i ) then
          call dswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
        end if
     end do
  else if ( incx > 1 ) then
     do i = k1, k2
        ip = ipiv( ix )
        if ( ip /= i ) then
          call dswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
        end if
        ix = ix + incx
     end do
  else if ( incx < 0 ) then
     do 30 i = k2, k1, -1
        ip = ipiv( ix )
        if ( ip /= i ) then
          call dswap( n, a( i, 1 ), lda, a( ip, 1 ), lda )
        end if
        ix = ix + incx
   30    continue
  end if

  return
end
function dnrm2 ( n, x, incx )
!
!***********************************************************************
!
!! DNRM2 returns the euclidean norm of a vector.
!
!
  real ( kind = 8 )      ONE         , zero
  PARAMETER           ( ONE = 1.0D+0, zero = 0.0D+0 )
!
  integer ( kind = 4 )                           INCX, N
  real ( kind = 8 )                  X( * )
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 )               IX
  real ( kind = 8 )      ABSXI, NORM, SCALE, SSQ
!
  if ( N < 1 .OR. INCX.LT.1 ) then
     NORM  = zero
  else if ( N == 1 ) then
     NORM  = ABS( X( 1 ) )
  else
     SCALE = zero
     SSQ   = ONE

     DO IX = 1, 1 + ( N - 1 )*INCX, INCX
        if ( X( IX ) /= zero ) then
           ABSXI = ABS( X( IX ) )
           if ( SCALE < ABSXI ) then
              SSQ   = ONE   + SSQ*( SCALE/ABSXI )**2
              SCALE = ABSXI
           else
              SSQ   = SSQ   +     ( ABSXI/SCALE )**2
           end if
        end if
     end do
     NORM  = SCALE * SQRT( SSQ )
  end if

  DNRM2 = NORM
  return
end
subroutine dorg2r( m, n, k, a, lda, tau, work, info )
!
!***********************************************************************
!
!! DORG2R ???
!
  integer ( kind = 4 )            INFO, K, LDA, M, N
  real ( kind = 8 )   A( LDA, * ), TAU( * ), WORK( * )
!
!  DORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) integer
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) real ( kind = 8 ) array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.
!
!  LDA     (input) integer
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) real ( kind = 8 ) array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace) real ( kind = 8 ) array, dimension (N)
!
!  INFO    (output) integer
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value
!
!
  real ( kind = 8 )   ONE, zero
  PARAMETER          ( ONE = 1.0D+0, zero = 0.0D+0 )
!     ..
!     .. Local Scalars ..
  integer ( kind = 4 )            I, J, L
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLARF, DSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
!
!     Test the input arguments
!
  INFO = 0
  if ( M < 0 ) then
     INFO = -1
  else if ( N < 0 .OR. N > M ) then
     INFO = -2
  else if ( K < 0 .OR. K > N ) then
     INFO = -3
  else if ( LDA < MAX( 1, M ) ) then
     INFO = -5
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DORG2R', -INFO )
     return
  end if
!
!     Quick return if possible
!
  if ( N <= 0 ) then
    return
  end if
!
!     Initialise columns k+1:n to columns of the unit matrix
!
  DO 20 J = K + 1, N
     DO L = 1, M
        A( L, J ) = zero
     end do
     A( J, J ) = ONE
   20 continue
!
  DO 40 I = K, 1, -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
     if ( I < N ) then
        A( I, I ) = ONE
        CALL DLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                       A( I, I+1 ), LDA, WORK )
     end if
     if ( I < M ) then
       CALL DSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
     end if
     A( I, I ) = ONE - TAU( I )
!
!        Set A(1:i-1,i) to zero
!
     DO 30 L = 1, I - 1
        A( L, I ) = zero
   30    continue
   40 continue
  return
end
subroutine dorgqr( m, n, k, a, lda, tau, work, lwork, info )
!
!***********************************************************************
!
!! DORGQR ???
!
  integer ( kind = 4 )            INFO, K, LDA, LWORK, M, N
  real ( kind = 8 )   A( LDA, * ), TAU( * ), WORK( LWORK )
!
!  DORGQR generates an M-by-N real matrix Q with orthonormal columns,
!  which is defined as the first N columns of a product of K elementary
!  reflectors of order M
!
!        Q  =  H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF.
!
!  Arguments
!
!  M       (input) integer
!          The number of rows of the matrix Q. M >= 0.
!
!  N       (input) integer
!          The number of columns of the matrix Q. M >= N >= 0.
!
!  K       (input) integer
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.
!
!  A       (input/output) real ( kind = 8 ) array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the M-by-N matrix Q.
!
!  LDA     (input) integer
!          The first dimension of the array A. LDA >= max(1,M).
!
!  TAU     (input) real ( kind = 8 ) array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  WORK    (workspace/output) real ( kind = 8 ) array, dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) integer
!          The dimension of the array WORK. LWORK >= max(1,N).
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!  INFO    (output) integer
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument has an illegal value
!
  real ( kind = 8 )   zero
  PARAMETER          ( zero = 0.0D+0 )
!
  integer ( kind = 4 )            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK, NB
  integer ( kind = 4 ) NBMIN, NX
!     ..
!     .. External Subroutines ..
  EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
  integer ( kind = 4 )            ILAENV
  EXTERNAL           ILAENV
!
!     Test the input arguments
!
  INFO = 0
  if ( M < 0 ) then
     INFO = -1
  else if ( N < 0 .OR. N > M ) then
     INFO = -2
  else if ( K < 0 .OR. K > N ) then
     INFO = -3
  else if ( LDA < MAX( 1, M ) ) then
     INFO = -5
  else if ( LWORK < MAX( 1, N ) ) then
     INFO = -8
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DORGQR', -INFO )
     return
  end if
!
!     Quick return if possible
!
  if ( N <= 0 ) then
     WORK( 1 ) = 1
     return
  end if
!
!     Determine the block size.
!
  NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
  NBMIN = 2
  NX = 0
  IWS = N
  if ( NB > 1 .and. NB < K ) then
!
!        Determine when to cross over from blocked to unblocked code.
!
     NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
     if ( NX < K ) then
!
!           Determine if workspace is large enough for blocked code.
!
        LDWORK = N
        IWS = LDWORK*NB
        if ( LWORK < IWS ) then
!
!              Not enough workspace to use optimal NB:  reduce NB and
!              determine the minimum value of NB.
!
           NB = LWORK / LDWORK
           NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
        end if
     end if
  end if
!
  if ( NB >= NBMIN .and. NB < K .AND. NX.LT.K ) then
!
!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.
!
     KI = ( ( K-NX-1 ) / NB )*NB
     KK = MIN( K, KI+NB )
!
!        Set A(1:kk,kk+1:n) to zero.
!
     DO 20 J = KK + 1, N
        DO 10 I = 1, KK
           A( I, J ) = zero
   10       continue
   20    continue
  else
     KK = 0
  end if
!
!     Use unblocked code for the last or only block.
!
  if ( KK < N ) then
    CALL DORG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA, &
                   TAU( KK+1 ), WORK, IINFO )
  end if

  if ( KK > 0 ) then
!
!        Use blocked code
!
     DO 50 I = KI + 1, 1, -NB
        IB = MIN( NB, K-I+1 )
        if ( I+IB <= N ) then
!
!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)
!
           CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB, &
                           A( I, I ), LDA, TAU( I ), WORK, LDWORK )
!
!              Apply H to A(i:m,i+ib:n) from the left
!
           CALL DLARFB( 'Left', 'No transpose', 'Forward', &
                           'Columnwise', M-I+1, N-I-IB+1, IB, &
                           A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ), &
                           LDA, WORK( IB+1 ), LDWORK )
        end if
!
!  Apply H to rows i:m of current block
!
        CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK, IINFO )
!
!  Set rows 1:i-1 of current block to zero
!
        DO 40 J = I, I + IB - 1
           DO L = 1, I - 1
              A( L, J ) = zero
           end do
   40       continue
   50    continue
  end if

  WORK( 1 ) = IWS
  return
end
subroutine dscal ( n, da, dx, incx )
!
!*******************************************************************************
!
!! DSCAL scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx  <=  0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
  real ( kind = 8 ) da,dx(*)
  integer ( kind = 4 ) i,incx,m,n,nincx
!
  if (  n <= 0 .or. incx <= 0 )return
  if ( incx == 1)go to 20
!
!        code for increment not equal to 1
!
  nincx = n*incx
  do i = 1,nincx,incx
    dx(i) = da*dx(i)
  end do
  return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 continue

  m = mod(n,5)
  if (  m  ==  0 ) go to 40
  dx(1:m) = da*dx(1:m)
  if (  n  <  5 ) return

40 continue

  do i = m+1,n,5
    dx(i) = da*dx(i)
    dx(i + 1) = da*dx(i + 1)
    dx(i + 2) = da*dx(i + 2)
    dx(i + 3) = da*dx(i + 3)
    dx(i + 4) = da*dx(i + 4)
  end do

  return
end
subroutine dswap ( n, x, incx, y, incy )
!
!*******************************************************************************
!
!! DSWAP interchanges two vectors.
!
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    08 April 1999
!
!  Reference:
!
!    Lawson, Hanson, Kincaid, Krogh,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    Algorithm 539,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, pages 308-323.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vectors.
!
!    Input/output, real ( kind = 8 ) X(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCX, the increment between successive entries of X.
!
!    Input/output, real ( kind = 8 ) Y(*), one of the vectors to swap.
!
!    Input, integer ( kind = 4 ) INCY, the increment between successive elements of Y.
!
  integer ( kind = 4 ) i
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) incy
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  real ( kind = 8 ) temp
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) y(*)
!
  if ( n <= 0 ) then

  else if ( incx == 1 .and. incy == 1 ) then

    m = mod ( n, 3 )

    do i = 1, m
      temp = x(i)
      x(i) = y(i)
      y(i) = temp
    end do

    do i = m+1, n, 3

      temp = x(i)
      x(i) = y(i)
      y(i) = temp

      temp = x(i+1)
      x(i+1) = y(i+1)
      y(i+1) = temp

      temp = x(i+2)
      x(i+2) = y(i+2)
      y(i+2) = temp

    end do

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    if ( incy >= 0 ) then
      iy = 1
    else
      iy = ( - n + 1 ) * incy + 1
    end if

    do i = 1, n
      temp = x(ix)
      x(ix) = y(iy)
      y(iy) = temp
      ix = ix + incx
      iy = iy + incy
    end do

  end if

  return
end
subroutine dtbsv(uplo,trans,diag,n,k,a,lda,x,incx)
!
!***********************************************************************
!
!! DTBSV ???
!
  integer ( kind = 4 ) incx, k, lda, n
  character        diag, trans, uplo
  real ( kind = 8 )   a( lda, * ), x( * )
!
!  DTBSV  solves one of the systems of equations
!
!     A*x = b,   or   A'*x = b,
!
!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular band matrix, with ( k + 1 )
!  diagonals.
!
!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.
!
!  Parameters
!
!  UPLO   - character*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:
!
!              TRANS = 'N' or 'n'   A*x = b.
!
!              TRANS = 'T' or 't'   A'*x = b.
!
!              TRANS = 'C' or 'c'   A'*x = b.
!
!           Unchanged on exit.
!
!  DIAG   - character*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  K      - integer.
!           On entry with UPLO = 'U' or 'u', K specifies the number of
!           super-diagonals of the matrix A.
!           On entry with UPLO = 'L' or 'l', K specifies the number of
!           sub-diagonals of the matrix A.
!           K must satisfy  0  <=  K.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!           by n part of the array A must contain the upper triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row
!           ( k + 1 ) of the array, the first super-diagonal starting at
!           position 2 in row k, and so on. The top left k by k triangle
!           of the array A is not referenced.
!           The following program segment will transfer an upper
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = K + 1 - J
!                    DO 10, I = MAX( 1, J - K ), J
!                       A( M + I, J ) = matrix( I, J )
!              10    continue
!              20 continue
!
!           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!           by n part of the array A must contain the lower triangular
!           band part of the matrix of coefficients, supplied column by
!           column, with the leading diagonal of the matrix in row 1 of
!           the array, the first sub-diagonal starting at position 1 in
!           row 2, and so on. The bottom right k by k triangle of the
!           array A is not referenced.
!           The following program segment will transfer a lower
!           triangular band matrix from conventional full matrix storage
!           to band storage:
!
!                 DO 20, J = 1, N
!                    M = 1 - J
!                    DO 10, I = J, MIN( N, J + K )
!                       A( M + I, J ) = matrix( I, J )
!              10    continue
!              20 continue
!
!           Note that when DIAG = 'U' or 'u' the elements of the array A
!           corresponding to the diagonal elements of the matrix are not
!           referenced, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           ( k + 1 ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.
!
!  INCX   - integer.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!     .. Parameters ..
  real ( kind = 8 )   zero
  parameter        ( zero = 0.0d+0 )
!     .. Local Scalars ..
  real ( kind = 8 )   temp
  integer ( kind = 4 )            i, info, ix, j, jx, kplus1, kx, l
  logical            nounit
  logical            lsame
!
!     Test the input parameters.
!
  info = 0
  if     ( .not.lsame( uplo , 'u' ).and. &
              .not.lsame( uplo , 'l' )      ) then
     info = 1
  else if ( .not.lsame( trans, 'n' ).and. &
              .not.lsame( trans, 't' ).and. &
              .not.lsame( trans, 'c' )      ) then
     info = 2
  else if ( .not.lsame( diag , 'u' ).and. &
              .not.lsame( diag , 'n' )      ) then
     info = 3
  else if ( n < 0 ) then
     info = 4
  else if ( k < 0 ) then
     info = 5
  else if ( lda < ( k + 1 ) ) then
     info = 7
  else if ( incx == 0 ) then
     info = 9
  end if
  if ( info /= 0 ) then
     call xerbla( 'dtbsv ', info )
     return
  end if
!
!     Quick return if possible.
!
  if ( n == 0 ) then
    return
  end if

  nounit = lsame( diag, 'n' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
  if ( incx <= 0 ) then
     kx = 1 - ( n - 1 )*incx
  else if ( incx /= 1 ) then
     kx = 1
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
  if ( lsame( trans, 'n' ) ) then
!
!        Form  x : =  inv( A )*x.
!
     if ( lsame( uplo, 'u' ) ) then
        kplus1 = k + 1
        if ( incx == 1 ) then
           do 20, j = n, 1, -1
              if ( x( j ) /= zero ) then
                 l = kplus1 - j
                 if ( nounit ) then
                   x( j ) = x( j )/a( kplus1, j )
                 end if
                 temp = x( j )
                 do 10, i = j - 1, max( 1, j - k ), -1
                    x( i ) = x( i ) - temp*a( l + i, j )
   10                continue
              end if
   20          continue
        else
           kx = kx + ( n - 1 )*incx
           jx = kx
           do 40, j = n, 1, -1
              kx = kx - incx
              if ( x( jx ) /= zero ) then
                 ix = kx
                 l  = kplus1 - j
                 if ( nounit ) then
                   x( jx ) = x( jx )/a( kplus1, j )
                 end if
                 temp = x( jx )
                 do 30, i = j - 1, max( 1, j - k ), -1
                    x( ix ) = x( ix ) - temp*a( l + i, j )
                    ix      = ix      - incx
   30                continue
              end if
              jx = jx - incx
   40          continue
        end if
     else
        if ( incx == 1 ) then
           do 60, j = 1, n
              if ( x( j ) /= zero ) then
                 l = 1 - j
                 if ( nounit ) then
                   x( j ) = x( j )/a( 1, j )
                 end if
                 temp = x( j )
                 do 50, i = j + 1, min( n, j + k )
                    x( i ) = x( i ) - temp*a( l + i, j )
   50                continue
              end if
   60          continue
        else
           jx = kx
           do 80, j = 1, n
              kx = kx + incx
              if ( x( jx ) /= zero ) then
                 ix = kx
                 l  = 1  - j
                 if ( nounit ) then
                   x( jx ) = x( jx )/a( 1, j )
                 end if
                 temp = x( jx )
                 do 70, i = j + 1, min( n, j + k )
                    x( ix ) = x( ix ) - temp*a( l + i, j )
                    ix      = ix      + incx
   70                continue
              end if
              jx = jx + incx
   80          continue
        end if
     end if
  else
!
!        Form  x : =  inv( A')*x.
!
     if ( lsame( uplo, 'u' ) ) then
        kplus1 = k + 1
        if ( incx == 1 ) then
           do 100, j = 1, n
              temp = x( j )
              l    = kplus1 - j
              do 90, i = max( 1, j - k ), j - 1
                 temp = temp - a( l + i, j )*x( i )
   90             continue
              if ( nounit ) then
                temp = temp/a( kplus1, j )
              end if
              x( j ) = temp
  100          continue
        else
           jx = kx
           do 120, j = 1, n
              temp = x( jx )
              ix   = kx
              l    = kplus1  - j
              do 110, i = max( 1, j - k ), j - 1
                 temp = temp - a( l + i, j )*x( ix )
                 ix   = ix   + incx
  110             continue
              if ( nounit ) then
                temp = temp/a( kplus1, j )
              end if
              x( jx ) = temp
              jx      = jx   + incx
              if ( j > k ) then
                kx = kx + incx
              end if
  120          continue
        end if
     else
        if ( incx == 1 ) then
           do 140, j = n, 1, -1
              temp = x( j )
              l    = 1      - j
              do 130, i = min( n, j + k ), j + 1, -1
                 temp = temp - a( l + i, j )*x( i )
  130             continue
              if ( nounit ) then
                temp = temp/a( 1, j )
              end if
              x( j ) = temp
  140          continue
        else
           kx = kx + ( n - 1 )*incx
           jx = kx
           do 160, j = n, 1, -1
              temp = x( jx )
              ix   = kx
              l    = 1       - j
              do 150, i = min( n, j + k ), j + 1, -1
                 temp = temp - a( l + i, j )*x( ix )
                 ix   = ix   - incx
  150             continue
              if ( nounit ) then
                temp = temp/a( 1, j )
              end if
              x( jx ) = temp
              jx      = jx   - incx
              if ( ( n - j ) >= k ) then
                kx = kx - incx
              end if
  160          continue
        end if
     end if
  end if

  return
end
subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb )
!
!***********************************************************************
!
!! DTRMM ???
!
!
  character       SIDE, UPLO, TRANSA, DIAG
  integer ( kind = 4 )            M, N, LDA, LDB
  real ( kind = 8 )   ALPHA
  real ( kind = 8 )   A( LDA, * ), B( LDB, * )
!
!  DTRMM  performs one of the matrix-matrix operations
!
!     B : =  alpha*op( A )*B,   or   B := alpha*B*op( A ),
!
!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  Parameters
!
!  SIDE   - character*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:
!
!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!
!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
!
!           Unchanged on exit.
!
!  UPLO   - character*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - character*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - character*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - integer.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - real ( kind = 8 ) array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.
!
!  LDB    - integer.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
  logical            LSAME
  EXTERNAL           LSAME
!
  EXTERNAL           XERBLA

!
  logical            LSIDE, NOUNIT, UPPER
  integer ( kind = 4 )            I, INFO, J, K, NROWA
  real ( kind = 8 )   TEMP
!
  real ( kind = 8 )   ONE         , zero
  PARAMETER        ( ONE = 1.0D+0, zero = 0.0D+0 )
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
  LSIDE  = LSAME( SIDE  , 'L' )
  if ( LSIDE ) then
     NROWA = M
  else
     NROWA = N
  end if
  NOUNIT = LSAME( DIAG  , 'N' )
  UPPER  = LSAME( UPLO  , 'U' )
!
  INFO   = 0
  if (      ( .NOT.LSIDE                ).and. &
              ( .NOT.LSAME( SIDE  , 'R' ) )      ) then
     INFO = 1
  else if ( ( .NOT.UPPER                ).and. &
             ( .NOT.LSAME( UPLO  , 'L' ) )      ) then
     INFO = 2
  else if ( ( .NOT.LSAME( TRANSA, 'N' ) ).and. &
              ( .NOT.LSAME( TRANSA, 'T' ) ).and. &
              ( .NOT.LSAME( TRANSA, 'C' ) )      ) then
     INFO = 3
  else if ( ( .NOT.LSAME( DIAG  , 'U' ) ).and. &
              ( .NOT.LSAME( DIAG  , 'N' ) )      ) then
     INFO = 4
  else if ( M   < 0               ) then
     INFO = 5
  else if ( N   < 0               ) then
     INFO = 6
  else if ( LDA < MAX( 1, NROWA ) ) then
     INFO = 9
  else if ( LDB < MAX( 1, M     ) ) then
     INFO = 11
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DTRMM ', INFO )
     return
  end if
!
!     Quick return if possible.
!
  if ( N == 0 ) then
    return
  end if
!
!     And when  alpha == zero.
!
  if ( ALPHA == zero ) then
     DO 20, J = 1, N
        DO 10, I = 1, M
           B( I, J ) = zero
   10       continue
   20    continue
     return
  end if
!
!     Start the operations.
!
  if ( LSIDE ) then
     if ( LSAME( TRANSA, 'N' ) ) then
!
!           Form  B : =  alpha*A*B.
!
        if ( UPPER ) then
           DO 50, J = 1, N
              DO 40, K = 1, M
                 if ( B( K, J ).NE.zero ) then
                    TEMP = ALPHA*B( K, J )
                    DO I = 1, K - 1
                       B( I, J ) = B( I, J ) + TEMP*A( I, K )
                    end do
                    if ( NOUNIT ) then
                      TEMP = TEMP*A( K, K )
                    end if
                    B( K, J ) = TEMP
                 end if
   40             continue
   50          continue
        else
           DO 80, J = 1, N
              DO 70 K = M, 1, -1
                 if ( B( K, J ).NE.zero ) then
                    TEMP      = ALPHA*B( K, J )
                    B( K, J ) = TEMP
                    if ( NOUNIT ) then
                      B( K, J ) = B( K, J )*A( K, K )
                    end if
                    DO 60, I = K + 1, M
                       B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   continue
                 end if
   70             continue
   80          continue
        end if
     else
!
!           Form  B : =  alpha*B*A'.
!
        if ( UPPER ) then
           DO 110, J = 1, N
              DO 100, I = M, 1, -1
                 TEMP = B( I, J )
                 if ( NOUNIT ) then
                   TEMP = TEMP*A( I, I )
                 end if
                 DO 90, K = 1, I - 1
                    TEMP = TEMP + A( K, I )*B( K, J )
   90                continue
                 B( I, J ) = ALPHA*TEMP
  100             continue
  110          continue
        else
           DO 140, J = 1, N
              DO 130, I = 1, M
                 TEMP = B( I, J )
                 if ( NOUNIT ) then
                   TEMP = TEMP*A( I, I )
                 end if
                 DO 120, K = I + 1, M
                    TEMP = TEMP + A( K, I )*B( K, J )
  120                continue
                 B( I, J ) = ALPHA*TEMP
  130             continue
  140          continue
        end if
     end if
  else
     if ( LSAME( TRANSA, 'N' ) ) then
!
!           Form  B : =  alpha*B*A.
!
        if ( UPPER ) then
           DO 180, J = N, 1, -1
              TEMP = ALPHA
              if ( NOUNIT ) then
                TEMP = TEMP*A( J, J )
              end if
              DO 150, I = 1, M
                 B( I, J ) = TEMP*B( I, J )
  150             continue
              DO 170, K = 1, J - 1
                 if ( A( K, J ).NE.zero ) then
                    TEMP = ALPHA*A( K, J )
                    DO 160, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   continue
                 end if
  170             continue
  180          continue
        else
           DO 220, J = 1, N
              TEMP = ALPHA
              if ( NOUNIT ) then
                TEMP = TEMP*A( J, J )
              end if
              DO 190, I = 1, M
                 B( I, J ) = TEMP*B( I, J )
  190             continue
              DO 210, K = J + 1, N
                 if ( A( K, J ).NE.zero ) then
                    TEMP = ALPHA*A( K, J )
                    DO 200, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   continue
                 end if
  210             continue
  220          continue
        end if
     else
!
!           Form  B : =  alpha*B*A'.
!
        if ( UPPER ) then
           DO 260, K = 1, N
              DO 240, J = 1, K - 1
                 if ( A( J, K ).NE.zero ) then
                    TEMP = ALPHA*A( J, K )
                    DO 230, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   continue
                 end if
  240             continue
              TEMP = ALPHA
              if ( NOUNIT ) then
                TEMP = TEMP*A( K, K )
              end if
              if ( TEMP.NE.ONE ) then
                 DO 250, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  250                continue
              end if
  260          continue
        else
           DO 300, K = N, 1, -1
              DO 280, J = K + 1, N
                 if ( A( J, K ).NE.zero ) then
                    TEMP = ALPHA*A( J, K )
                    DO 270, I = 1, M
                       B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   continue
                 end if
  280             continue
              TEMP = ALPHA
              if ( NOUNIT ) then
                TEMP = TEMP*A( K, K )
              end if
              if ( TEMP.NE.ONE ) then
                 DO 290, I = 1, M
                    B( I, K ) = TEMP*B( I, K )
  290                continue
              end if
  300          continue
        end if
     end if
  end if

  return
end
subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx )
!
!***********************************************************************
!
!! DTRMV ???
!
!
  integer ( kind = 4 )            INCX, LDA, N
  character       DIAG, TRANS, UPLO
  real ( kind = 8 )   A( LDA, * ), X( * )
!
!  DTRMV  performs one of the matrix-vector operations
!
!     x : =  A*x,   or   x := A'*x,
!
!  where x is an n element vector and  A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.
!
!  Parameters
!
!  UPLO   - character*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANS  - character*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:
!
!              TRANS = 'N' or 'n'   x := A*x.
!
!              TRANS = 'T' or 't'   x := A'*x.
!
!              TRANS = 'C' or 'c'   x := A'*x.
!
!           Unchanged on exit.
!
!  DIAG   - character*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.
!
!  X      - real ( kind = 8 ) array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.
!
!  INCX   - integer.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
  real ( kind = 8 )   zero
  PARAMETER        ( zero = 0.0D+0 )
!     .. Local Scalars ..
  real ( kind = 8 )   TEMP
  integer ( kind = 4 )            I, INFO, IX, J, JX, KX
  logical            NOUNIT
!     .. External Functions ..
  logical            LSAME
  EXTERNAL           LSAME
!     .. External Subroutines ..
  EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
!
!     Test the input parameters.
!
  INFO = 0
  IF     ( .NOT.LSAME( UPLO , 'U' ).and. &
             .NOT.LSAME( UPLO , 'L' )      ) then
     INFO = 1
  else if ( .NOT.LSAME( TRANS, 'N' ).and. &
              .NOT.LSAME( TRANS, 'T' ).and. &
              .NOT.LSAME( TRANS, 'C' )      ) then
     INFO = 2
  else if ( .NOT.LSAME( DIAG , 'U' ).and. &
              .NOT.LSAME( DIAG , 'N' )      ) then
     INFO = 3
  else if ( N < 0 ) then
     INFO = 4
  else if ( LDA < MAX( 1, N ) ) then
     INFO = 6
  else if ( INCX == 0 ) then
     INFO = 8
  end if
  if ( INFO.NE.0D+00 ) then
     CALL XERBLA( 'DTRMV ', INFO )
     return
  end if

  if ( N == 0 ) then
    return
  end if

  NOUNIT = LSAME( DIAG, 'N' )
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
  if ( INCX <= 0 ) then
     KX = 1 - ( N - 1 )*INCX
  else if ( INCX.NE.1 ) then
     KX = 1
  end if
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
  if ( LSAME( TRANS, 'N' ) ) then
!
!        Form  x : =  A*x.
!
     if ( LSAME( UPLO, 'U' ) ) then
        if ( INCX == 1 ) then
           DO 20, J = 1, N
              if ( X( J ).NE.zero ) then
                 TEMP = X( J )
                 DO 10, I = 1, J - 1
                    X( I ) = X( I ) + TEMP*A( I, J )
   10                continue
                 if ( NOUNIT ) then
                   X( J ) = X( J )*A( J, J )
                 end if
              end if
   20          continue
        else
           JX = KX
           DO 40, J = 1, N
              if ( X( JX ).NE.zero ) then
                 TEMP = X( JX )
                 IX   = KX
                 DO 30, I = 1, J - 1
                    X( IX ) = X( IX ) + TEMP*A( I, J )
                    IX      = IX      + INCX
   30                continue
                 if ( NOUNIT ) then
                   X( JX ) = X( JX )*A( J, J )
                 end if
              end if
              JX = JX + INCX
   40          continue
        end if
     else
        if ( INCX == 1 ) then
           DO 60, J = N, 1, -1
              if ( X( J ).NE.zero ) then
                 TEMP = X( J )
                 DO 50, I = N, J + 1, -1
                    X( I ) = X( I ) + TEMP*A( I, J )
   50                continue
                 if ( NOUNIT ) then
                   X( J ) = X( J )*A( J, J )
                 end if
              end if
   60          continue
        else
           KX = KX + ( N - 1 )*INCX
           JX = KX
           DO 80, J = N, 1, -1
              if ( X( JX ).NE.zero ) then
                 TEMP = X( JX )
                 IX   = KX
                 DO 70, I = N, J + 1, -1
                    X( IX ) = X( IX ) + TEMP*A( I, J )
                    IX      = IX      - INCX
   70                continue
                 if ( NOUNIT ) then
                   X( JX ) = X( JX )*A( J, J )
                 end if
              end if
              JX = JX - INCX
   80          continue
        end if
     end if
  else
!
!        Form  x : =  A'*x.
!
     if ( LSAME( UPLO, 'U' ) ) then
        if ( INCX == 1 ) then
           DO 100, J = N, 1, -1
              TEMP = X( J )
              if ( NOUNIT ) then
                TEMP = TEMP*A( J, J )
              end if
              DO 90, I = J - 1, 1, -1
                 TEMP = TEMP + A( I, J )*X( I )
   90             continue
              X( J ) = TEMP
  100          continue
        else
           JX = KX + ( N - 1 )*INCX
           DO 120, J = N, 1, -1
              TEMP = X( JX )
              IX   = JX
              if ( NOUNIT ) then
                TEMP = TEMP*A( J, J )
              end if
              DO 110, I = J - 1, 1, -1
                 IX   = IX   - INCX
                 TEMP = TEMP + A( I, J )*X( IX )
  110             continue
              X( JX ) = TEMP
              JX      = JX   - INCX
  120          continue
        end if
     else
        if ( INCX == 1 ) then
           DO J = 1, N
              TEMP = X( J )
              if ( NOUNIT ) then
                TEMP = TEMP*A( J, J )
              end if
              DO I = J + 1, N
                 TEMP = TEMP + A( I, J )*X( I )
              end do
              X( J ) = TEMP
           end do
        else
           JX = KX
           DO 160, J = 1, N
              TEMP = X( JX )
              IX   = JX
              if ( NOUNIT ) then
                 TEMP = TEMP*A( J, J )
              end if
              DO 150, I = J + 1, N
                 IX   = IX   + INCX
                 TEMP = TEMP + A( I, J )*X( IX )
  150             continue
              X( JX ) = TEMP
              JX      = JX   + INCX
  160          continue
        end if
     end if
  end if

  return
end
subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
!
!***********************************************************************
!
!! DTRSM ???
!
  character        side, uplo, transa, diag
  integer ( kind = 4 )            m, n, lda, ldb
  real ( kind = 8 )   alpha
  real ( kind = 8 )   a( lda, * ), b( ldb, * )
!
!  DTRSM  solves one of the matrix equations
!
!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!
!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!
!     op( A ) = A   or   op( A ) = A'.
!
!  The matrix X is overwritten on B.
!
!  Parameters
!
!  SIDE   - character*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:
!
!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!
!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!
!           Unchanged on exit.
!
!  UPLO   - character*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:
!
!              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!
!              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!
!           Unchanged on exit.
!
!  TRANSA - character*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:
!
!              TRANSA = 'N' or 'n'   op( A ) = A.
!
!              TRANSA = 'T' or 't'   op( A ) = A'.
!
!              TRANSA = 'C' or 'c'   op( A ) = A'.
!
!           Unchanged on exit.
!
!  DIAG   - character*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:
!
!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!
!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.
!
!           Unchanged on exit.
!
!  M      - integer.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.
!
!  N      - integer.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.
!
!  ALPHA  - real ( kind = 8 ).
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.
!
!  A      - real ( kind = 8 ) array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.
!
!  LDA    - integer.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.
!
!  B      - real ( kind = 8 ) array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.
!
!  LDB    - integer.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.
!
  logical            lsame
  external           lsame
  external           xerbla
  logical            lside, nounit, upper
  integer ( kind = 4 )            i, info, j, k, nrowa
  real ( kind = 8 )   temp
!     .. Parameters ..
  real ( kind = 8 )   one         , zero
  parameter        ( one = 1.0d+0, zero = 0.0d+0 )
!
!     Test the input parameters.
!
  lside  = lsame( side  , 'l' )
  if ( lside ) then
     nrowa = m
  else
     nrowa = n
  end if
  nounit = lsame( diag  , 'n' )
  upper  = lsame( uplo  , 'u' )
!
  info   = 0
  if (      ( .not.lside                ).and. &
             ( .not.lsame( side  , 'r' ) )      ) then
     info = 1
  else if ( ( .not.upper                ).and. &
    ( .not.lsame( uplo  , 'l' ) )      ) then
     info = 2
  else if ( ( .not.lsame( transa, 'n' ) ).and. &
              ( .not.lsame( transa, 't' ) ).and. &
              ( .not.lsame( transa, 'c' ) )      ) then
     info = 3
  else if ( ( .not.lsame( diag  , 'u' ) ).and. &
              ( .not.lsame( diag  , 'n' ) )      ) then
     info = 4
  else if ( m   < 0               ) then
     info = 5
  else if ( n   < 0               ) then
     info = 6
  else if ( lda < max( 1, nrowa ) ) then
     info = 9
  else if ( ldb < max( 1, m     ) ) then
     info = 11
  end if
  if ( info /= 0 ) then
     call xerbla( 'dtrsm ', info )
     return
  end if

  if ( n == 0 ) then
    return
  end if
!
!     And when  alpha == zero.
!
  if ( alpha == zero ) then
     do 20, j = 1, n
        do 10, i = 1, m
           b( i, j ) = zero
   10       continue
   20    continue
     return
  end if
!
!     Start the operations.
!
  if ( lside ) then
     if ( lsame( transa, 'n' ) ) then
!
!           Form  B : =  alpha*inv( A )*B.
!
        if ( upper ) then
           do 60, j = 1, n
              if ( alpha /= one ) then
                 do 30, i = 1, m
                    b( i, j ) = alpha*b( i, j )
   30                continue
              end if
              do 50, k = m, 1, -1
                 if ( b( k, j ) /= zero ) then
                    if ( nounit ) then
                          b( k, j ) = b( k, j )/a( k, k )
                    end if
                    do 40, i = 1, k - 1
                       b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   40                   continue
                 end if
   50             continue
   60          continue
        else
           do 100, j = 1, n
              if ( alpha /= one ) then
                 do i = 1, m
                    b( i, j ) = alpha*b( i, j )
                 end do
              end if
              do 90 k = 1, m
                 if ( b( k, j ) /= zero ) then
                    if ( nounit ) then
                      b( k, j ) = b( k, j )/a( k, k )
                    end if
                    do i = k + 1, m
                       b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
                    end do
                 end if
   90             continue
  100          continue
        end if
     else
!
!           Form  B : =  alpha*inv( A' )*B.
!
        if ( upper ) then
           do 130, j = 1, n
              do 120, i = 1, m
                 temp = alpha*b( i, j )
                 do 110, k = 1, i - 1
                    temp = temp - a( k, i )*b( k, j )
  110                continue
                 if ( nounit ) then
                   temp = temp/a( i, i )
                 end if
                 b( i, j ) = temp
  120             continue
  130          continue
        else
           do 160, j = 1, n
              do 150, i = m, 1, -1
                 temp = alpha*b( i, j )
                 do 140, k = i + 1, m
                    temp = temp - a( k, i )*b( k, j )
  140                continue
                 if ( nounit ) then
                   temp = temp/a( i, i )
                 end if
                 b( i, j ) = temp
  150             continue
  160          continue
        end if
     end if
  else
     if ( lsame( transa, 'n' ) ) then
!
!           Form  B : =  alpha*B*inv( A ).
!
        if ( upper ) then
           do 210, j = 1, n
              if ( alpha /= one ) then
                 do 170, i = 1, m
                    b( i, j ) = alpha*b( i, j )
  170                continue
              end if
              do 190, k = 1, j - 1
                 if ( a( k, j ) /= zero ) then
                    do 180, i = 1, m
                       b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  180                   continue
                 end if
  190             continue
              if ( nounit ) then
                 temp = one/a( j, j )
                 do i = 1, m
                    b( i, j ) = temp*b( i, j )
                 end do
              end if
  210          continue
        else
           do 260, j = n, 1, -1
              if ( alpha /= one ) then
                 do i = 1, m
                    b( i, j ) = alpha*b( i, j )
                 end do
              end if
              do 240, k = j + 1, n
                 if ( a( k, j ) /= zero ) then
                    do 230, i = 1, m
                       b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  230                   continue
                 end if
  240             continue
              if ( nounit ) then
                 temp = one/a( j, j )
                 do 250, i = 1, m
                   b( i, j ) = temp*b( i, j )
  250                continue
              end if
  260          continue
        end if
     else
!
!           Form  B : =  alpha*B*inv( A' ).
!
        if ( upper ) then
           do 310, k = n, 1, -1
              if ( nounit ) then
                 temp = one/a( k, k )
                 do 270, i = 1, m
                    b( i, k ) = temp*b( i, k )
  270                continue
              end if
              do 290, j = 1, k - 1
                 if ( a( j, k ) /= zero ) then
                    temp = a( j, k )
                    do 280, i = 1, m
                       b( i, j ) = b( i, j ) - temp*b( i, k )
  280                   continue
                 end if
  290             continue
              if ( alpha /= one ) then
                 do 300, i = 1, m
                    b( i, k ) = alpha*b( i, k )
  300                continue
              end if
  310          continue
        else
           do 360, k = 1, n
              if ( nounit ) then
                 temp = one/a( k, k )
                 do 320, i = 1, m
                    b( i, k ) = temp*b( i, k )
  320                continue
              end if
              do 340, j = k + 1, n
                 if ( a( j, k ) /= zero ) then
                    temp = a( j, k )
                    do 330, i = 1, m
                       b( i, j ) = b( i, j ) - temp*b( i, k )
  330                   continue
                 end if
  340             continue
              if ( alpha /= one ) then
                 do 350, i = 1, m
                    b( i, k ) = alpha*b( i, k )
  350                continue
              end if
  360          continue
        end if
     end if
  end if
!
  return
end
function idamax(n,dx,incx)
!
!***********************************************************************
!
!! IDAMAX FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE.
!
  real ( kind = 8 ) dmax
  real ( kind = 8 ) dx(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idamax
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) ix
  integer ( kind = 4 ) n
!
  idamax = 0

  if ( n  <  1 ) return

  idamax = 1
  if ( n == 1)return

  if ( incx == 1)go to 20
!
  ix = 1
  dmax = abs(dx(1))
  ix = ix + incx

  do i = 2,n

    if ( abs(dx(ix)) > dmax) then
      idamax = i
      dmax = abs(dx(ix))
      ix = ix + incx
    end if

  end do

  return
!
   20 dmax = abs(dx(1))

  do i = 2,n
    if ( abs(dx(i)) > dmax) then
      idamax = i
      dmax = abs(dx(i))
    end if
  end do

  return
end
function ilaenv( ispec, name, opts, n1, n2, n3, n4 )
!
!***********************************************************************
!
!! ILAENV ???
!
!
  integer ( kind = 4 ) ilaenv
  character*( * )    NAME, OPTS
  integer ( kind = 4 )            ISPEC, N1, N2, N3, N4
!
!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.
!
!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.
!
!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.
!
!  Arguments
!
!  ISPEC   (input) integer
!          Specifies the parameter to be returned as the value of
!          ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.
!
!  NAME    (input) character*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) character*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) integer
!  N2      (input) integer
!  N3      (input) integer
!  N4      (input) integer
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) integer
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!
!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:
!
!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      if ( NB <= 1 ) NB = MAX( 1, N )
!
!     .. Local Scalars ..
  logical            CNAME, SNAME
  character       C1
  character ( len = 2 )        C2, C4
  character ( len = 3 )       C3
  character ( len = 6 )       SUBNAM
  integer ( kind = 4 )            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
!
  go to ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
  ILAENV = -1
  return
!
  100 continue
!
!     Convert NAME to upper case if the first character is lower case.
!
  ILAENV = 1
  SUBNAM = NAME
  IC = ICHAR( SUBNAM( 1:1 ) )
  IZ = ICHAR( 'Z' )
  if ( IZ == 90 .OR. IZ == 122 ) then
!
!        ASCII character set
!
     if ( IC >= 97 .and. IC <= 122 ) then
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO 10 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           if ( IC >= 97 .and. IC <= 122 ) then
                 SUBNAM( I:I ) = CHAR( IC-32 )
           end if
   10       continue
     end if
!
  else if ( IZ == 233 .OR. IZ == 169 ) then
!
!        EBCDIC character set
!
     if ( ( IC >= 129 .and. IC <= 137 ) .OR. &
            ( IC >= 145 .and. IC <= 153 ) .OR. &
            ( IC >= 162 .and. IC <= 169 ) ) then
        SUBNAM( 1:1 ) = CHAR( IC+64 )
        DO 20 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           if ( ( IC >= 129 .and. IC <= 137 ) .OR. &
                  ( IC >= 145 .and. IC <= 153 ) .OR. &
                  ( IC >= 162 .and. IC <= 169 ) ) then
                 SUBNAM( I:I ) = CHAR( IC+64 )
           end if
   20       continue
     end if
!
  else if ( IZ == 218 .OR. IZ == 250 ) then
!
!        Prime machines:  ASCII+128
!
     if ( IC >= 225 .and. IC <= 250 ) then
        SUBNAM( 1:1 ) = CHAR( IC-32 )
        DO 30 I = 2, 6
           IC = ICHAR( SUBNAM( I:I ) )
           if ( IC >= 225 .and. IC <= 250 ) then
                 SUBNAM( I:I ) = CHAR( IC-32 )
           end if
   30       continue
     end if
  end if
!
  C1 = SUBNAM( 1:1 )
  SNAME = C1 == 'S' .OR. C1 == 'D'
  CNAME = C1 == 'C' .OR. C1 == 'Z'
  if ( .NOT.( CNAME .OR. SNAME ) ) then
    return
  end if
  C2 = SUBNAM( 2:3 )
  C3 = SUBNAM( 4:6 )
  C4 = C3( 2:3 )
!
  go to ( 110, 200, 300 ) ISPEC
!
  110 continue
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or real ( kind = 8 ).
!
  NB = 1
!
  if ( C2 == 'GE' ) then
     if ( C3 == 'TRF' ) then
        if ( SNAME ) then
           NB = 64
        else
           NB = 64
        end if
     else if ( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. &
                 C3 == 'QLF' ) then
        if ( SNAME ) then
           NB = 32
        else
           NB = 32
        end if
     else if ( C3 == 'HRD' ) then
        if ( SNAME ) then
           NB = 32
        else
           NB = 32
        end if
     else if ( C3 == 'BRD' ) then
        if ( SNAME ) then
           NB = 32
        else
           NB = 32
        end if
     else if ( C3 == 'TRI' ) then
        if ( SNAME ) then
           NB = 64
        else
           NB = 64
        end if
     end if
  else if ( C2 == 'PO' ) then
     if ( C3 == 'TRF' ) then
        if ( SNAME ) then
           NB = 64
        else
           NB = 64
        end if
     end if
  else if ( C2 == 'SY' ) then
     if ( C3 == 'TRF' ) then
        if ( SNAME ) then
           NB = 64
        else
           NB = 64
        end if
     else if ( SNAME .and. C3 == 'TRD' ) then
        NB = 1
     else if ( SNAME .and. C3 == 'GST' ) then
        NB = 64
     end if
  else if ( CNAME .and. C2 == 'HE' ) then
     if ( C3 == 'TRF' ) then
        NB = 64
     else if ( C3 == 'TRD' ) then
        NB = 1
     else if ( C3 == 'GST' ) then
        NB = 64
     end if
  else if ( SNAME .and. C2 == 'OR' ) then
     if ( C3( 1:1 ) == 'G' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NB = 32
        end if
     else if ( C3( 1:1 ) == 'M' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NB = 32
        end if
     end if
  else if ( CNAME .and. C2 == 'UN' ) then
     if ( C3( 1:1 ) == 'G' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NB = 32
        end if
     else if ( C3( 1:1 ) == 'M' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NB = 32
        end if
     end if
  else if ( C2 == 'GB' ) then
     if ( C3 == 'TRF' ) then
        if ( SNAME ) then
           if ( N4 <= 64 ) then
              NB = 1
           else
              NB = 32
           end if
        else
           if ( N4 <= 64 ) then
              NB = 1
           else
              NB = 32
           end if
        end if
     end if
  else if ( C2 == 'PB' ) then
     if ( C3 == 'TRF' ) then
        if ( SNAME ) then
           if ( N2 <= 64 ) then
              NB = 1
           else
              NB = 32
           end if
        else
           if ( N2 <= 64 ) then
              NB = 1
           else
              NB = 32
           end if
        end if
     end if
  else if ( C2 == 'TR' ) then
     if ( C3 == 'TRI' ) then
        if ( SNAME ) then
           NB = 64
        else
           NB = 64
        end if
     end if
  else if ( C2 == 'LA' ) then
     if ( C3 == 'UUM' ) then
        if ( SNAME ) then
           NB = 64
        else
           NB = 64
        end if
     end if
  else if ( SNAME .and. C2 == 'ST' ) then
     if ( C3 == 'EBZ' ) then
        NB = 1
     end if
  end if
  ILAENV = NB
  return
!
  200 continue
!
!     ISPEC = 2:  minimum block size
!
  NBMIN = 2
  if ( C2 == 'GE' ) then
     if ( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. C3 == 'QLF' ) then
        if ( SNAME ) then
           NBMIN = 2
        else
           NBMIN = 2
        end if
     else if ( C3 == 'HRD' ) then
        if ( SNAME ) then
           NBMIN = 2
        else
           NBMIN = 2
        end if
     else if ( C3 == 'BRD' ) then
        if ( SNAME ) then
           NBMIN = 2
        else
           NBMIN = 2
        end if
     else if ( C3 == 'TRI' ) then
        if ( SNAME ) then
           NBMIN = 2
        else
           NBMIN = 2
        end if
     end if
  else if ( C2 == 'SY' ) then
     if ( C3 == 'TRF' ) then
        if ( SNAME ) then
           NBMIN = 8
        else
           NBMIN = 8
        end if
     else if ( SNAME .and. C3 == 'TRD' ) then
        NBMIN = 2
     end if
  else if ( CNAME .and. C2 == 'HE' ) then
     if ( C3 == 'TRD' ) then
        NBMIN = 2
     end if
  else if ( SNAME .and. C2 == 'OR' ) then
     if ( C3( 1:1 ) == 'G' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NBMIN = 2
        end if
     else if ( C3( 1:1 ) == 'M' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NBMIN = 2
        end if
     end if
  else if ( CNAME .and. C2 == 'UN' ) then
     if ( C3( 1:1 ) == 'G' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NBMIN = 2
        end if
     else if ( C3( 1:1 ) == 'M' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
             C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
             C4 == 'BR' ) then
           NBMIN = 2
        end if
     end if
  end if
  ILAENV = NBMIN
  return
!
  300 continue
!
!     ISPEC = 3:  crossover point
!
  NX = 0
  if ( C2 == 'GE' ) then
     if ( C3 == 'QRF' .OR. C3 == 'RQF' .OR. C3 == 'LQF' .OR. C3 == 'QLF' ) then
        if ( SNAME ) then
           NX = 128
        else
           NX = 128
        end if
     else if ( C3 == 'HRD' ) then
        if ( SNAME ) then
           NX = 128
        else
           NX = 128
        end if
     else if ( C3 == 'BRD' ) then
        if ( SNAME ) then
           NX = 128
        else
           NX = 128
        end if
     end if
  else if ( C2 == 'SY' ) then
     if ( SNAME .and. C3 == 'TRD' ) then
        NX = 1
     end if
  else if ( CNAME .and. C2 == 'HE' ) then
     if ( C3 == 'TRD' ) then
        NX = 1
     end if
  else if ( SNAME .and. C2 == 'OR' ) then
     if ( C3( 1:1 ) == 'G' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
               C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
               C4 == 'BR' ) then
           NX = 128
        end if
     end if
  else if ( CNAME .and. C2 == 'UN' ) then
     if ( C3( 1:1 ) == 'G' ) then
        if ( C4 == 'QR' .OR. C4 == 'RQ' .OR. C4 == 'LQ' .OR. &
               C4 == 'QL' .OR. C4 == 'HR' .OR. C4 == 'TR' .OR. &
               C4 == 'BR' ) then
           NX = 128
        end if
     end if
  end if
  ILAENV = NX
  return
!
  400 continue
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
  ILAENV = 6
  return
!
  500 continue
!
!     ISPEC = 5:  minimum column dimension (not used)
!
  ILAENV = 2
  return
!
  600 continue
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
  ILAENV = INT( real( MIN( N1, N2 ) )*1.6E0 )
  return
!
  700 continue
!
!     ISPEC = 7:  number of processors (not used)
!
  ILAENV = 1
  return
!
  800 continue
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
  ILAENV = 50
  return
end
function lsame( ca, cb )
!
!***********************************************************************
!
!! LSAME ???
!
  character          CA, CB
  logical lsame
!     ..
!
!  Purpose
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!
!  CA      (input) character*1
!  CB      (input) character*1
!          CA and CB specify the single characters to be compared.
!
  integer ( kind = 4 )            INTA, INTB, ZCODE
!
!     Test if the characters are equal
!
  LSAME = CA == CB
  if ( LSAME ) then
    return
  end if
!
!     Now test for equivalence if both characters are alphabetic.
!
  ZCODE = ICHAR( 'Z' )
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
  INTA = ICHAR( CA )
  INTB = ICHAR( CB )
!
  if ( ZCODE == 90 .OR. ZCODE == 122 ) then
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
     if ( INTA >= 97 .and. INTA <= 122 ) INTA = INTA - 32
     if ( INTB >= 97 .and. INTB <= 122 ) INTB = INTB - 32
!
  else if ( ZCODE == 233 .OR. ZCODE == 169 ) then
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
     if ( INTA >= 129 .and. INTA <= 137 .OR. &
            INTA >= 145 .and. INTA <= 153 .OR. &
            INTA >= 162 .and. INTA <= 169 ) INTA = INTA + 64
     if ( INTB >= 129 .and. INTB <= 137 .OR. &
             INTB >= 145 .and. INTB <= 153 .OR. &
             INTB >= 162 .and. INTB <= 169 ) INTB = INTB + 64
!
  else if ( ZCODE == 218 .OR. ZCODE == 250 ) then
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
     if ( INTA >= 225 .and. INTA <= 250 ) INTA = INTA - 32
     if ( INTB >= 225 .and. INTB <= 250 ) INTB = INTB - 32
  end if
  LSAME = INTA == INTB

  return
end
subroutine xerbla(srname,info)
!
!***********************************************************************
!
!! XERBLA is an error handler for the LAPACK routines.
!
!
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!
!  SRNAME  (input) character*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) integer
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
  character ( len = 6 )       SRNAME
  integer ( kind = 4 )            INFO
!
  WRITE( *, FMT = 9999 )SRNAME, INFO
  STOP

 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
   'an illegal value' )
end
function dlamch ( CMACH )
!
!***********************************************************************
!
!! DLAMCH ???
!
  character          CMACH
!
!  dlamch determines real ( kind = 8 ) machine parameters.
!
!  Arguments
!
!  CMACH   (input) character
!          Specifies the value to be returned by dlamch:
!          = 'E' or 'e',   dlamch := eps
!          = 'S' or 's ,   dlamch := sfmin
!          = 'B' or 'b',   dlamch := base
!          = 'P' or 'p',   dlamch := eps*base
!          = 'N' or 'n',   dlamch := t
!          = 'R' or 'r',   dlamch := rnd
!          = 'M' or 'm',   dlamch := emin
!          = 'U' or 'u',   dlamch := rmin
!          = 'L' or 'l',   dlamch := emax
!          = 'O' or 'o',   dlamch := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0D+00 when rounding occurs in addition, 0.0D+00 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
  real ( kind = 8 )   ONE, ZERO
  PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!
  real ( kind = 8 ) dlamch
  logical            FIRST, LRND
  integer ( kind = 4 )            BETA, IMAX, IMIN, IT
  real ( kind = 8 )   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN
  real ( kind = 8 )       RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
  logical            LSAME
  EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
  EXTERNAL           dlamc2
!     ..
!     .. Save statement ..
  SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN
  save EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
  DATA               FIRST / .TRUE. /
!
  if ( FIRST ) THEN
     FIRST = .FALSE.
     CALL dlamc2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
     BASE = BETA
     T = IT
     if ( LRND ) THEN
        RND = ONE
        EPS = ( BASE**( 1-IT ) ) / 2
     else
        RND = ZERO
        EPS = BASE**( 1-IT )
     END IF
     PREC = EPS*BASE
     EMIN = IMIN
     EMAX = IMAX
     SFMIN = RMIN
     SMALL = ONE / RMAX
     if ( SMALL >= SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
        SFMIN = SMALL*( ONE+EPS )
     END IF
  END IF
!
  if ( LSAME( CMACH, 'E' ) ) THEN
     RMACH = EPS
  else if ( LSAME( CMACH, 'S' ) ) THEN
     RMACH = SFMIN
  else if ( LSAME( CMACH, 'B' ) ) THEN
     RMACH = BASE
  else if ( LSAME( CMACH, 'P' ) ) THEN
     RMACH = PREC
  else if ( LSAME( CMACH, 'N' ) ) THEN
     RMACH = T
  else if ( LSAME( CMACH, 'R' ) ) THEN
     RMACH = RND
  else if ( LSAME( CMACH, 'M' ) ) THEN
     RMACH = EMIN
  else if ( LSAME( CMACH, 'U' ) ) THEN
     RMACH = RMIN
  else if ( LSAME( CMACH, 'L' ) ) THEN
     RMACH = EMAX
  else if ( LSAME( CMACH, 'O' ) ) THEN
     RMACH = RMAX
  END IF

  dlamch = RMACH
  RETURN
end
subroutine dlamc1( BETA, T, RND, IEEE1 )
!
!***********************************************************************
!
!! DLAMC1 ???
!
  logical            IEEE1, RND
  integer ( kind = 4 )            BETA, T
!
!  dlamc1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!
!  BETA    (output) integer
!          The base of the machine.
!
!  T       (output) integer
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) logical
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) logical
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!
!  The routine is based on the routine  ENVRON  by Malcolm and
!  incorporates suggestions by Gentleman and Marovich. See
!
!     Malcolm M. A. (1972) Algorithms to reveal properties of
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!        that reveal properties of floating point arithmetic units.
!        Comms. of the ACM, 17, 276-277.
!
!     .. Local Scalars ..
  logical            FIRST, LIEEE1, LRND
  integer ( kind = 4 )            LBETA, LT
  real ( kind = 8 )   A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
  real ( kind = 8 )   dlamc3
  EXTERNAL           dlamc3
!     ..
!     .. Save statement ..
  SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
  DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
  if ( FIRST ) THEN
     FIRST = .FALSE.
     ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  dlamc3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0D+00 **m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0D+00 ) = a.
!
     A = 1
     C = 1

  do while ( C == ONE )
        A = 2*A
        C = dlamc3( A, ONE )
        C = dlamc3( C, -A )
  end do
!
!        Now compute  b = 2.0D+00 **m  with the smallest positive integer m
!        such that
!
!           fl( a + b )  >  a.
!
     B = 1
     C = dlamc3( A, B )

  do while ( c == a )
        B = 2*B
        C = dlamc3( A, B )
  end do
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
     QTR = ONE / 4
     SAVEC = C
     C = dlamc3( C, -A )
     LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
     B = LBETA
     F = dlamc3( B / 2, -B / 100 )
     C = dlamc3( F, A )
     if ( C == A ) THEN
        LRND = .TRUE.
     else
        LRND = .FALSE.
     END IF
     F = dlamc3( B / 2, B / 100 )
     C = dlamc3( F, A )
     if ( ( LRND ) .AND. ( C == A ) ) then
       LRND = .FALSE.
     end if
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
     T1 = dlamc3( B / 2, A )
     T2 = dlamc3( B / 2, SAVEC )
     LIEEE1 = ( T1 == A ) .AND. ( T2 > SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer ( kind = 4 ) part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0D+00 ) = 1.0.
!
     LT = 0
     A = 1
     C = 1

     do while ( C == ONE )
        LT = LT + 1
        A = A*LBETA
        C = dlamc3( A, ONE )
        C = dlamc3( C, -A )
     end do

  END IF

  BETA = LBETA
  T = LT
  RND = LRND
  IEEE1 = LIEEE1
  RETURN
end
subroutine dlamc2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!***********************************************************************
!
!! DLAMC2 ???
!
  logical            RND
  integer ( kind = 4 )            BETA, EMAX, EMIN, T
  real ( kind = 8 )   EPS, RMAX, RMIN
!
!  dlamc2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!
!  BETA    (output) integer
!          The base of the machine.
!
!  T       (output) integer
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) logical
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) real ( kind = 8 )
!          The smallest positive number such that
!
!             fl( 1.0D+00 - EPS )  <  1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) integer
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) real ( kind = 8 )
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) integer
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) real ( kind = 8 )
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
  logical            FIRST, IEEE, IWARN, LIEEE1, LRND
  integer ( kind = 4 )            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT
  integer ( kind = 4 )               NGNMIN, NGPMIN
  real ( kind = 8 )   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE
  real ( kind = 8 )                   SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
  real ( kind = 8 )   dlamc3
  EXTERNAL           dlamc3
!     ..
!     .. External Subroutines ..
  EXTERNAL           dlamc1, DLAMC4, DLAMC5
!
  SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX
  save LRMIN, LT
!     ..
!     .. Data statements ..
  DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
  if ( FIRST ) THEN
     FIRST = .FALSE.
     ZERO = 0
     ONE = 1
     TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  dlamc3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        dlamc1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
     CALL dlamc1( LBETA, LT, LRND, LIEEE1 )
!
!        Start to find EPS.
!
     B = LBETA
     A = B**( -LT )
     LEPS = A
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
     B = TWO / 3
     HALF = ONE / 2
     SIXTH = dlamc3( B, -HALF )
     THIRD = dlamc3( SIXTH, SIXTH )
     B = dlamc3( THIRD, -HALF )
     B = dlamc3( B, SIXTH )
     B = ABS( B )
     if ( B < LEPS ) then
       B = LEPS
     end if

     LEPS = 1

     do WHILE( ( LEPS > B ).AND.( B > ZERO ) )
        LEPS = B
        C = dlamc3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
        C = dlamc3( HALF, -C )
        B = dlamc3( HALF, C )
        C = dlamc3( HALF, -B )
        B = dlamc3( HALF, C )
     end do

     if ( A < LEPS ) then
       LEPS = A
     end if
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
     RBASE = ONE / LBETA
     SMALL = ONE
     DO I = 1, 3
        SMALL = dlamc3( SMALL*RBASE, ZERO )
     end do
     A = dlamc3( ONE, SMALL )
     CALL dlamc4( NGPMIN, ONE, LBETA )
     CALL dlamc4( NGNMIN, -ONE, LBETA )
     CALL dlamc4( GPMIN, A, LBETA )
     CALL dlamc4( GNMIN, -A, LBETA )
     IEEE = .FALSE.
!
     if ( ( NGPMIN == NGNMIN ) .AND. ( GPMIN == GNMIN ) ) THEN
        if ( NGPMIN == GPMIN ) THEN
           LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
        else if ( ( GPMIN-NGPMIN ) == 3 ) THEN
           LEMIN = NGPMIN - 1 + LT
           IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
        else
           LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
           IWARN = .TRUE.
        END IF
!
     else if ( ( NGPMIN == GPMIN ) .AND. ( NGNMIN == GNMIN ) ) THEN
        if ( ABS( NGPMIN-NGNMIN ) == 1 ) THEN
           LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
        else
           LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
           IWARN = .TRUE.
        END IF

     else if ( ( ABS( NGPMIN-NGNMIN ) == 1 ) .AND.( GPMIN == GNMIN ) ) THEN
        if ( ( GPMIN-MIN( NGPMIN, NGNMIN ) ) == 3 ) THEN
           LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
        else
           LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
           IWARN = .TRUE.
        END IF
!
     else
        LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
        IWARN = .TRUE.
     END IF
!**
! Comment out this if block if EMIN is ok
     if ( IWARN ) THEN
        FIRST = .TRUE.
        WRITE( 6, FMT = 9999 )LEMIN
     END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine dlamc1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
     IEEE = IEEE .OR. LIEEE1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
     LRMIN = 1
     DO 30 I = 1, 1 - LEMIN
        LRMIN = dlamc3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call dlamc5 to compute EMAX and RMAX.
!
     CALL dlamc5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
  END IF
!
  BETA = LBETA
  T = LT
  RND = LRND
  EPS = LEPS
  EMIN = LEMIN
  RMIN = LRMIN
  EMAX = LEMAX
  RMAX = LRMAX
!
  RETURN

 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
           '  EMIN = ', I8, / &
           ' If, after inspection, the value EMIN looks', &
           ' acceptable please comment out ', &
           / ' the IF block as marked within the code of routine', &
           ' dlamc2,', / ' otherwise supply EMIN explicitly.', / )
end
function dlamc3( A, B )
!
!***********************************************************************
!
!! DLAMC3 is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!
!  A, B    (input) real ( kind = 8 )
!          The values A and B.
!
  real ( kind = 8 )   A, B
  real ( kind = 8 ) dlamc3
!
  dlamc3 = A + B

  RETURN
end
subroutine dlamc4( EMIN, START, BASE )
!
!***********************************************************************
!
!! DLAMC4 ???
!
  integer ( kind = 4 )            BASE, EMIN
  real ( kind = 8 )   START
!
!  dlamc4 is a service routine for DLAMC2.
!
!  Arguments
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) real ( kind = 8 )
!          The starting point for determining EMIN.
!
!  BASE    (input) integer
!          The base of the machine.
!
  integer ( kind = 4 )            I
  real ( kind = 8 )   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
  real ( kind = 8 )   dlamc3
  EXTERNAL           dlamc3
!     ..
!     .. Executable Statements ..
!
  A = START
  ONE = 1
  RBASE = ONE / BASE
  ZERO = 0
  EMIN = 1
  B1 = dlamc3( A*RBASE, ZERO )
  C1 = A
  C2 = A
  D1 = A
  D2 = A
   10 CONTINUE

  if ( ( C1 == A ) .AND. ( C2 == A ) .AND. ( D1 == A ) .AND.( D2 == A ) ) THEN
     EMIN = EMIN - 1
     A = B1
     B1 = dlamc3( A / BASE, ZERO )
     C1 = dlamc3( B1*BASE, ZERO )
     D1 = ZERO
     DO I = 1, BASE
        D1 = D1 + B1
     end do
     B2 = dlamc3( A*RBASE, ZERO )
     C2 = dlamc3( B2 / RBASE, ZERO )
     D2 = ZERO
     DO I = 1, BASE
        D2 = D2 + B2
     end do
     GO TO 10
  END IF
!+    END WHILE
!
  RETURN
end
subroutine dlamc5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!*****************************************************************************80
!
!! DLAMC5 ???
!

  logical            IEEE
  integer ( kind = 4 )            BETA, EMAX, EMIN, P
  real ( kind = 8 )   RMAX
!
!  dlamc5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!
!  BETA    (input) integer
!          The base of floating-point arithmetic.
!
!  P       (input) integer
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) integer
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) logical
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) integer
!          The largest exponent before overflow
!
!  RMAX    (output) real ( kind = 8 )
!          The largest machine floating-point number.
!
!     .. Parameters ..
  real ( kind = 8 )   ZERO, ONE
  PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
!     ..
!     .. Local Scalars ..
  integer ( kind = 4 )            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
  real ( kind = 8 )   OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
  real ( kind = 8 )   dlamc3
  EXTERNAL           dlamc3
!     ..
!     .. Intrinsic Functions ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
  LEXP = 1
  EXBITS = 1
   10 CONTINUE
  TRY = LEXP*2
  if ( TRY <= ( -EMIN ) ) THEN
     LEXP = TRY
     EXBITS = EXBITS + 1
     GO TO 10
  END IF
  if ( LEXP == -EMIN ) THEN
     UEXP = LEXP
  else
     UEXP = TRY
     EXBITS = EXBITS + 1
  END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
  if ( ( UEXP+EMIN ) > ( -LEXP-EMIN ) ) THEN
     EXPSUM = 2*LEXP
  else
     EXPSUM = 2*UEXP
  END IF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
  EMAX = EXPSUM + EMIN - 1
  NBITS = 1 + EXBITS + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
  if ( ( MOD( NBITS, 2 ) == 1 ) .AND. ( BETA == 2 ) ) THEN
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
     EMAX = EMAX - 1
  END IF
!
  if ( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
     EMAX = EMAX - 1
  END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0D+00 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0D+00 - BETA**(-P), being careful that the
!     result is less than 1.0D+00 .
!
  RECBAS = ONE / BETA
  Z = BETA - ONE
  Y = ZERO
  DO 20 I = 1, P
     Z = Z*RECBAS
     if ( Y < ONE ) then
       OLDY = Y
     end if
     Y = dlamc3( Y, Z )
   20 CONTINUE

  if ( Y >= ONE ) then
    Y = OLDY
  end if
!
!     Now multiply by BETA**EMAX to get RMAX.
!
  DO 30 I = 1, EMAX
     Y = dlamc3( Y*BETA, ZERO )
   30 CONTINUE

  RMAX = Y
  RETURN
end
subroutine dlarfg(n,alpha,x,incx,tau)
!
!*****************************************************************************80
!
!! DLARFG ???
!
!
!  DLARFG generates a real elementary reflector H of order n, such
!  that
!
!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )
!
!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form
!
!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )
!
!  where tau is a real scalar and v is a real (n-1)-element
!  vector.
!
!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.
!
!  Otherwise  1 <= tau <= 2.
!
!  Arguments
!
!  N       (input) integer
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) real ( kind = 8 )
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) real ( kind = 8 ) array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) integer
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) real ( kind = 8 )
!          The value tau.
!
!
  real ( kind = 8 ) one
  real ( kind = 8 ) zero

  parameter (one = 1.0d0)
  parameter (zero = 0.0d0)
!
  real ( kind = 8 ) alpha
  real ( kind = 8 ) beta
  real ( kind = 8 ) dlamch
  real ( kind = 8 ) dlapy2
  real ( kind = 8 ) dnrm2
  integer ( kind = 4 ) incx
  integer ( kind = 4 ) j
  integer ( kind = 4 ) knt
  integer ( kind = 4 ) n
  real ( kind = 8 ) rsafmn
  real ( kind = 8 ) safmin
  real ( kind = 8 ) tau
  real ( kind = 8 ) x(*)
  real ( kind = 8 ) xnorm
!
  if ( n <= 1) then
     TAU = zero
     return
  end if
!
  XNORM = DNRM2( N-1, X, INCX )
!
  if ( XNORM == zero ) then
!
!        H  =  I
!
     TAU = zero
  else
!
!        general case
!
     BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
     SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )

!
!  XNORM, BETA may be inaccurate; scale X and recompute them
!
     if ( ABS( BETA ) < SAFMIN ) then
        RSAFMN = ONE / SAFMIN
        KNT = 0
   10       continue
        knt = knt + 1
        call DSCAL( N-1, RSAFMN, X, INCX )
        BETA = BETA*RSAFMN
        ALPHA = ALPHA*RSAFMN
        if ( ABS( BETA ) < SAFMIN ) go to 10
!
!  New BETA is at most 1, at least SAFMIN
!
        XNORM = DNRM2( N-1, X, INCX )
        BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
        TAU = ( BETA-ALPHA ) / BETA
        CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!  If ALPHA is subnormal, it may lose relative accuracy
!
        ALPHA = BETA
        DO J = 1, KNT
           ALPHA = ALPHA*SAFMIN
        end do
     else
        TAU = ( BETA-ALPHA ) / BETA
        CALL DSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
        ALPHA = BETA
     end if
  end if

  return
end
