program main

!*****************************************************************************80
!
!! MAIN is the main program for ALSCAL.
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     a l s c a l  pc  version                     last change: febr. 90
!
!
!     a l s c a l 8 4                                 version 84.1
!
!     ALSCAL performs metric or nonmetric:
!       multidimensional scaling (simple, replicated, weighted,
!         generalized, asymmetric, or weighted asymmetric models).
!       multidimensional unfolding (simple, replicated or weighted)
!
!     ALSCAL analyzes data which are either
!       two or three way
!       symmetric, asymmetric or rectangular
!       row or matrix conditional, or unconditional
!       replicated or unreplicated
!       binary, nominal, ordinal, interval or ratio
!       discrete or continuous
!       with or without missing data (any pattern)
!
!     The alternating least squares algorithm is described in:
!
!     1) Y Takane, F W Young, and  J Deleeuw, 
!        nonmetric individual differences multidimensional scaling. 
!        psychometrica,
!        vol. 42, no. 1, march 1977, pp. 7-67.
!
!     2) F W Young, Takane, y. and R Lewyckyj,
!        Three notes on ALSCAL,
!        psychometrika, vol. 43, 1978, p 433-435.
!
!     3) F W Young, 
!        an asymmetric Euclidean model for multiprocess asymmetric data.  
!        us-japan mds seminar proceedings, 1975.
!
!     4) F W Young,
!        Enhancements in alscal-82.  
!        sugi proceedings,
!        1982, 7, 633-642.
!
!     5) F W Young,
!        the general Euclidean model. 
!        in: law, h.g., et.al. (eds.), 
!        three-mode models for data analysis. 
!        praeger 1984
!
!
!     version 1.01           june 1974           Yoshio Takane
!     version 2.01           august 1974         Yoshio Takane
!     version 2.02           february 1976       Forrest W. Young
!     version 2.03f          october 1976        Forrest W. Young
!     version 2.03d          august 1976         Rostyslaw Lewyckyj
!     version 2.04           december 1976       Forrest W. Young
!     version 2.05           september 1977      Yoshio Takane
!     version 3.01f          december 1977       Rostyslaw Lewyckyj
!     version 4.01f          november 1978       Forrest W. Young
!     version 79.1sas        january 1979        Forrest W. Young
!     version 4.02f          february 1980       Forrest W. Young
!     version 79.5sas        december 1981       Forrest W. Young
!     version 82.1sas        may 1982            Forrest W. Young
!     version 83.1f          september 1982      jeff brooks
!     version 82.2sas        december 1982       Forrest W. Young
!     version 82.5sas        january  1983       Forrest W. Young
!     version 84.1f          june 1983           Forrest W. Young
!     version 84.2f/pc       february 1989       bernd erichson
!     version 84.3f/pc       february 1990       bernd erichson
!
!-pc--------------------------------------------------------------------
!     adapted for pc by
!       bernd erichson and alfred bischoff
!       herrngartenstr. 9, d-8501 kalchreuth, fed. rep. of germany
!-----------------------------------------------------------------------
!
!  this is the main routine for the fixed core version of alscal 84.1
!
!  this routine calls step0 which calculates the number of words
!  of storage required by the problem.  it then checks to see if
!  the program has sufficient memory allocated.  if not the program
!  terminates with a message indicating the amount of memory needed.
!  if sufficient memory is allocated this routine calls driver
!  which proceeds with the analysis of the data.
!
!-----------------------------------------------------------------------
!
!  program limits:
!
!    a) the number of dimensions of the multidimensional solution may
!       never exceed six.
!
!    b) all nominal or ordinal observations must be not greater than
!       9.0e20.
!
!    c) the maximum total problem size may not exceed the value set
!       for variable maxsiz below.
!
!-----------------------------------------------------------------------
!
!  i/o device assignments
!
!    all i/o devices are specified in the main routine.
!    the program uses 3 scratch files in addition to the 3 standard
!    units.   the standard units are assumed to be numbered 1, 2, and
!    3, and the scratch files 11, 12, 13.  these specifications
!    may be changed in the main routine.
!
!-----------------------------------------------------------------------
!
!  variables used in this routine
!
!     area     this array contains the work area for all dimensioned
!                  variables in the remaining routines.  the length
!                  of this array determines the maximum problem size.
!
!     maxsiz   the maximum problem size in words
!
!     prsize   the size of the problem being analyzed in words.
!
!     eoj    end of job indicator.
!
!     ionums   common area with input, output and scratch unit numbers.
!
!-----------------------------------------------------------------------
!
!  instructions for setting maximum problem size
!
!     a) choose an estimated maximum problem size, in words.
!        consultation with your computer center, plus knowledge of
!        the nature of your data, plus calculating the value of the
!        approximate problem size formula given in the users manual
!        will help in determing this value.
!        note that this value determine the amount of core, in words,
!        needed by alscal.
!
!     b) set the length of array area equal to your selected maximum
!        problem size.
!     c) set variable maxsiz equal to the value you have selected.
!     d) recompile the main routine, and relink alscal.
!     e) analyze your data.  if the maximum problem size is too small,
!        see below.

!-----------------------------------------------------------------------
!
!  instructions for changing maximum problem size
!
!     a) determine maximum problem size in words.  this is equal to the
!        estimated problem size printed by alscal when the maximum
!        problem size is too small to permit a problem to run.
!
!     b) change the length of array area to be at least as long as the
!        estimated problem size.  to be on the safe side, add
!        at least 10 percent to the size estimated by alscal.
!
!     c) change the value of maxsiz to equal the length of area.
!     d) recompile main and relink alscal.
!     e) you are now ready to try to reanalyze your data.
!
  integer, parameter :: maxsiz = 60000

  integer area(maxsiz)
  character ( len = 8 ) date
  integer eoj
  logical exist
  character ( len = 80 ) fname
  integer prsize
  character ( len = 10 ) time

  common /ionums/ in, nplt, lout, ndp, ndq, ndr, ndpp, indata

  call date_and_time ( date, time )
!
!  I/O unit specification
!
!  in=1 is the reader
!  nplt=2 is the punch
!  lout=3 is the printer
!  ndp=11 is a scratch unit
!  ndq=12 is a scratch unit
!  ndpp=13 is a scratch unit

  in = 1
  nplt = 2
  lout = 3
  ndp = 11
  ndq = 12
  ndpp = 13

  write ( *, * ) ' '
  write ( *, * ) 'ALSCALAL  AL        SCALALSC  ALSCALAL  ALSCALAL  AL      '
  write ( *, * ) 'LS    LS  LA        CA        LS        LS    LS  LA      '
  write ( *, * ) 'SC    SC  AL        AL        SC        SC    SC  AL      '
  write ( *, * ) 'CALALSCA  LS        LALSCALA  CA        CALALSCA  LS      '
  write ( *, * ) 'AL    AL  SC              AL  AL        AL    AL  SC      '
  write ( *, * ) 'LA    LA  CA              LS  LA        LA    LA  CA      '
  write ( *, * ) 'AL    AL  ALALSCAL  SCALALSC  ALSCALAL  AL    AL  ALALSCAL'
  write ( *, * ) ' '
  write ( *, * ) 'Alternating  least  squares   scaling    *    Version 84.1'
  write ( *, * ) 'by Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj'
  write ( *, * ) 'PC version by  B. Erichson and A. Bischoff  *  February 90'
  write ( *, * ) ' '
  write ( *, * ) 'Copyright 1977 by F W Young, Y Takane and R J Lewyckyj'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
  write ( *, * ) ' '
  write ( *, * ) 'Enter name of input file:'

  read ( *, '(a)' ) fname

  if ( len_trim ( fname ) == 0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'ALSCAL - Note.'
    write ( *, * ) '  Input file name was blank.'
    write ( *, * ) '  Normal end of execution.'
    stop
  end if

  exist = .false.

  inquire ( file = fname, exist = exist )

  if ( .not. exist ) then
    write ( *, * ) ' '
    write ( *, * ) 'ALSCAL - Fatal error!'
    write ( *, * ) '  Could not find the input file ', trim ( fname )
    write ( *, * ) '  End of execution.'
    stop
  end if

  open ( unit = in, file = fname, status = 'old' )

  write ( *, * ) ' '
  write ( *, * ) 'Enter the name of the output file:'

  read ( *, '(a)' ) fname

  if ( len_trim ( fname ) == 0 ) then
    fname = 'alscal.out'
    write ( *, * ) ' '
    write ( *, * ) 'The output file will be: ' // trim ( fname )
  end if

  open ( unit = lout, file = fname, status = 'replace' )

  open ( unit = ndp, status = 'scratch', form = 'unformatted' )

  open ( unit = ndq, status = 'scratch', form = 'unformatted' )

  open ( unit = ndpp, status = 'scratch', form = 'unformatted' )

  do

    eoj = - 1

    call step0 ( prsize, eoj )

    if ( eoj /= 0 ) then
      exit
    end if

    if ( prsize > maxsiz ) then
      write ( *, * ) ' '
      write ( *, * ) 'ALSCAL - Fatal error!'
      write ( *, * ) '  Not enough memory.'
      write ( *, * ) '  Recompile the main program after changing the dimension'
      write ( *, * ) '  of AREA to ', prsize
      write ( *, * ) '  and the value of MAXSIZ to ', prsize
      write ( *, * ) ' '
      write ( *, * ) '  Execution halted.'
      stop
    end if

    area(1:prsize) = 0

    write ( lout, * ) ' '
    write ( lout, * ) '  Maximum problem size :', maxsiz
    write ( lout, * ) '  Required problem size:', prsize

    call driver ( area, eoj )

    if ( eoj /= 0 ) then
      exit
    end if

  end do

  write ( *, * ) ' '
  write ( *, * ) 'ALSCAL - Note.'
  write ( *, * ) '  Normal end of execution.'

  stop
end
subroutine arnge ( w, x, ws, nb, ns, ndim, tr )

!*****************************************************************************80
!
!! ARNGE rearranges the weight and stimulus matrices.
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  integer nb
  integer ns
!
  integer j
  integer k
  integer nn(6)
  double precision t
  real tr(1)
  real w(ns,1)
  real ws(nb,1)
  real x(nb,1)
!
  common /block2/ ncst, nsim, nwe, ndmx, nab, ncol
!
  nn(1) = 1

  do j = 1, ndim

    t = 0.0d0

    if ( (nwe/2)*2 /= nwe ) then

      do i = 1, ns
        t = t + w(i,j)
      end do

    else if ( nwe >= 2 ) then

      do i = 1, nb
        t = t + ws(i,j)
      end do

    else

      do i = 1, nb
        t = t + x(i,j)**2
      end do

    end if

    tr(j) = - t

  end do
!
!  Sort TR into ascending order by insertion sort.
!
  do j = 2, ndim

    t = tr(j)
    k = j - 1

    do while ( k >= 1 )

      if ( t >= tr(k) ) then
        exit
      end if

      tr(k+1) = tr(k)
      nn(k+1) = nn(k)
      k = k - 1

    end do

    tr(k+1) = t
    nn(k+1) = j

  end do

  do j = 1, ndim

    k = nn(j)

    if ( j == k ) then
      cycle
    end if

    do i = 1, nb
      call r_swap ( x(i,j), x(i,k) )
    end do

    if ( (nwe/2)*2 /= nwe ) then

      do i = 1, ns
        call r_swap ( w(i,j), w(i,k) )
      end do

    end if

    if ( nwe >= 2 ) then

      do i = 1, nb
        call r_swap ( ws(i,j), ws(i,k) )
      end do

    end if

    do i = j+1, ndim

      if ( nn(i) == j ) then
        nn(i) = k
        exit
      end if

    end do

  end do

  return
end
subroutine bloc2 ( a, ll, iz, n )

!*****************************************************************************80
!
!! BLOC2 finds blocks of ties.
!
!  Modified:
!
!    09 May 2002
!
!  Author:
!
!    Forrest Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  integer n
!
  real a(n)
  integer iz
  integer ll(n)
!
  common /ionums/ in,nplt,lout,ndp,ndq,ndr,ndpp,indata
  common /block2/ ncst,nsim,nwe,ndmx,nab,ncol

  ii = 1

  do

    i = ii

    do

      ii = ii + 1

      if ( ii > n ) then
        exit
      end if

      if ( a(i) /= a(ii) ) then
        exit
      end if

    end do

    if ( ii - i /= 1 ) then

      if ( iz >= ndmx-1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ALSCAL - Fatal error!'
        write ( *, '(a,i6)' ) '  Maximum number of tie-blocks exceeds ', ndmx
        stop
      end if

      ll(iz) = i
      ll(iz+1) = ii - i
      iz = iz + 2

    end if

    if ( ii >= n ) then
      exit
    end if

  end do

  return
end
subroutine catses ( wa, wc, idy, n, civ, ncat )

!*****************************************************************************80
!
!! CATSES computes a discrete nominal transformation.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!        (disparities are means of all distances in category)
!
  real civ(1)
  integer idy(1)
  real wa(1)
  real wc(1)
!
  civ(1:ncat) = 0.0
  wc(1:ncat) = 0.0

!  At this moment wc is being used as scratch space.
!  Later it will be assigned new values to be returned to the
!  calling routine

  do i=1,n
    id=idy(i)
    civ(id)=civ(id)+wa(i)
    wc(id)=wc(id)+1.0
  end do

  do i=1,ncat
    civ(i) = civ(i) / wc(i)
  end do
!
!  now wc is assigned values to be returned to the caller
!
  do i = 1, n
    wc(i)=civ(idy(i))
  end do

  return
end
subroutine cjeig ( a, u, v, n, nd, b, w, alam, nft, nb )

!*****************************************************************************80
!
!! CJEIG computes the eigenvalues and eigenvectors of a real symmetric matrix.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     clint and jennings' method is used.
!        written by Yoshio Takane, 1974
!
  dimension a(n,1),u(n,1),v(n,1),b(nb,1),w(n,1),alam(1)
  double precision t1
!
  nd1=nd-1

  if ( nft == 1 ) then

    do j = 1, nd
      do i = 1, n
        u(i,j) = 0.0
      end do
      u(j,j) = 1.0
    end do

  end if

  do ll = 1, 30

    do i=1,n
      do j=1,nd
        t1=0.0
        do k=1,n
          t1=t1+a(i,k)*u(k,j)
        end do
        v(i,j)=t1
      end do
    end do

    do i=1,nd
      do j=1,nd
        t1=0.0
        do k=1,n
          t1=t1+u(k,i)*v(k,j)
        end do
        b(i,j)=t1
      end do
    end do

    call eigk ( b, alam, nd, nb )

    do 14 i=1,n
    do 14 j=1,nd
    t1=0.0
    do 41 k=1,nd
   41   t1=t1+v(i,k)*b(k,j)
   14   w(i,j)=t1


    do 15 i=1,nd
    do 15 j=1,nd

    t1=0.0
    do 51 k=1,n
   51   t1=t1+w(k,i)*w(k,j)

   15   b(i,j)=t1
!
!     cholesky factorization
!
    do 16 i=1,nd
    if ( b(i,i) <= 0.0) go to 18
    b(i,i)=sqrt(b(i,i))
    if ( i >= nd)go to 18
    j1=i+1
    do 16 j=j1,nd
    b(i,j)=b(i,j)/b(i,i)
    do 16 k=j1,j
   16   b(k,j)=b(k,j)-b(i,j)*b(i,k)
   18   continue
!
!     inverse of cholesky factor
!
    do i=1,nd
      b(i,i)=1.0/b(i,i)
      do j=1, i-1
        t1=0.0
        do k=j, i-1
          t1=t1-b(k,i)*b(k,j)
        end do
        b(i,j)=t1*b(i,i)
      end do
    end do

    v(1:n,1:nd) = w(1:n,1:nd)

    do i=1,n
      do j=1,nd
        t1=0.0
        do k=1,j
          t1=t1+v(i,k)*b(j,k)
        end do
        w(i,j)=t1
      end do
    end do
!
!  Test of convergence
!
    do 20 j=1,nd1
      do 20 i=1,n
        if ( abs(w(i,j)-u(i,j)) > 1.0e-5)go to 21
   20   continue

    u(1:n,1:nd) = w(1:n,1:nd)

    return

   21   continue

    u(1:n,1:nd) = w(1:n,1:nd)

  end do

  return
end
subroutine clear

!*****************************************************************************80
!
!! CLEAR clears the screen.
!
!
!  Use lahey system subroutine
!
!     call system ('cls')
!
!  Use assembler subroutine
!
!     call gclear
!
  return
end
subroutine coef ( u11, u12, u22, ub1, ub2, ndim, cfl, m, xn, nb, ndx )
!
!*****************************************************************************80
!
!! COEF obtains least stress coordinates.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     subroutine to obtain least stress coordinates
!     by solving a system of cubic equations for each coordinate
!     (Takane, Young and de leeuw, psychometrika, 1977)
!
  dimension u11(ndx,1),u12(ndx,1),u22(ndx,1),cfl(nb,1)
  dimension ub1(1),ub2(1),xn(1)
!
  do ll=1,30

    do i=1,ndim
      xn(i)=cfl(m,i)
    end do

    do l = 1, ndim

      a=2.0*u22(l,l)
      p=3.0*u12(l,l)/a
      q=u11(l,l)-2.0*ub2(l)
      r=-ub1(l)

      do i = 1, ndim
        if ( i /= l ) then
          q=q+2.0*cfl(m,i)*(u12(i,l)+u22(l,i)*cfl(m,i))
          r=r+cfl(m,i)*(u11(l,i)+u12(l,i)*cfl(m,i))
        end if
      end do

      q=q/a
      r=r/a
      call scube(p,q,r,cfl(m,l))

    end do

    do l=1,ndim
      if ( abs(cfl(m,l)-xn(l)) > 0.0001) then
        go to 20
      end if
    end do

    return

   20   continue

  end do

  return
end
subroutine distp ( w, cfl, x, wc, wd, ix, iy, iz, xx, ndsr, ws, nad )
!
!*****************************************************************************80
!
!! DISTP computes distances and disparities for coordinates and weights.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!  this routine computes distances and disparities (optimally scaled
!  data) for a given set of coordinates and weights by the Young,
!  de leeuw and Takane method
!
  double precision s
  integer ix(nt)
  integer iy(1)
  integer iz(1)
  dimension w(ns,1),cfl(nb,1),x(1),wc(1),wd(nt)
  dimension xx(1),ws(nb,1),ndsr(ns,nb),nad(ns,nb)
  dimension r(5,5),alph(5)

  common /block1/ nc,nd,big,nc2,ndt,nnc,nph,npt,nsc,epsi,ndim,ndx,ndxs, &
         ndxp,maxit,nadct,ndct,strso,strss,strss2,nb,ns,ndtyp,nps,&
         nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
!
!     compute squared distances
!
  last=0
  if ( nad(1,1) < 0)last=1
  if ( nad(1,1) < 0)nad(1,1)=-nad(1,1)
  osest=big
  strss=0.0
  lin=0
  if ( ndeg == 1.and.ncst == 1)lin=1
  n=0
  nx=2
  if ( nsim > 1) nx=1
  i1=nb
  do 37 l=1,ns
  do 35 i=nx,nb
  if ( nsim <= 1)i1=i-1
  do 35 j=1,i1
  n=n+1
  s=0.0
  if ( nwe-1)25,30,38
   25 do 26 k=1,ndim
   26 s=s+(cfl(i,k)-cfl(j,k))**2
  go to 34
   30 do 33 k=1,ndim
   33 s=s+w(l,k)*(cfl(i,k)-cfl(j,k))**2
  go to 34
   38 do 36 k=1,ndim
   36 s=s+w(l,k)*ws(i,k)*(cfl(i,k)-cfl(j,k))**2
   34 wc(n)=s
  wd(n)=s
   35 continue
   37 continue
!
!  Perform optimal scaling (compute disparities)
!
  n=1
  if ( ndtyp-3) 134,32,42

! -----  ordinal data  -----

   32 if ( nwc == 1) go to 434
  if ( nwc == 2) go to 506
!
!     unconditional ordinal data
!
  if ( nps /= 1) call prs(ix,iy,wc,nab,xx,iz,1)
  if ( nps == 1) call ses(ix,iy,wd,nab,iz,1)
  call trs(wd,x,ix,nab)
  if ( nab == nt)go to 150
  nab1=nab+1
  do 151 i=nab1,nt
  j=ix(i)
  if ( last == 1)osest=wc(j)
  151 x(j)=osest
  150 if ( last == 1)return
  call length(wc,x,nt)
  call mstrs(nt,x,wc,strss,strss2)
  strss=sqrt(strss)
  return
!
!     matrix conditional ordinal data
!
  434 do 435 l=1,ns
  nc3=nad(l,1)
  if ( nps /= 1) call prs(ix(n),iy,wc(n),nc3,xx,iz,l)
  if ( nps == 1) call ses(ix(n),iy,wd(n),nc3,iz,l)
  call trs(wd(n),x(n),ix(n),nc3)
  if ( nc3 >= nc2) go to 453
  ni=n+nc3
  ne=n+nc2-1
  do 152 i=ni,ne
  j=ix(i)+n-1
  if ( last == 1)osest=wc(j)
  152 x(j)=osest
  453 if ( last == 1)go to 1436
  call length(wc(n),x(n),nc2)
  call mstrs(nc2,x(n),wc(n),strss1,strss2)
  strss=strss+strss1
 1436 n=n+nc2
  435 continue
  strss=sqrt(strss/ns)
  return
!
!  Row conditional ordinal data
!
  506 n1=0
  do 507 l=1,ns
  do 507 i=1,nb
  nc3=nad(l,i)
  n1=n1+1
  if ( nsim > 3.and.i <= ncol)go to 504
  if ( nps /= 1) call prs(ix(n),iy,wc(n),nc3,xx,iz,n1)
  if ( nps == 1) call ses(ix(n),iy,wd(n),nc3,iz,n1)
  call trs(wd(n),x(n),ix(n),nc3)
  504 if ( nc3 >= nb)go to 505
  if ( nsim > 3.and.i <= ncol)x(n+i-1)=0.0
  ni=n+nc3
  ne=n+nb-1
  do 508 j=ni,ne
  k=ix(j)+n-1
  if ( last == 1)osest=wc(k)
  508 x(k)=osest
  505 if ( last == 1)go to 1506
  call length(wc(n),x(n),nb)
  call mstrs(nb,x(n),wc(n),strss1,strss2)
  strss=strss+strss1
 1506 n=n+nb
  507 continue
  na=nb
  if ( nsim > 3)na=nb-ncol
  strss=sqrt(strss/(na*ns))
  return
!
!  -----  nominal data  -----
!
   42 if ( nwc == 1) go to 436
  if ( nwc == 2) go to 509
!
!     unconditional nominal data
!
  call catses(wc,x,ix,nt,wd,ndct)
  if ( nab == nt.or.nps /= 1)go to 43
  do 201 i=1,nt
  if ( last == 1)osest=wc(i)
  201 if ( ix(i) == ndct) x(i)=osest
   43 if ( last == 1)return
  call length(wc,x,nt)
  call mstrs(nt,x,wc,strss,strss2)
  strss=sqrt(strss)
  return
!
!     matrix conditional nominal data
!
  436 do 437 l=1,ns
  call catses(wc(n),x(n),ix(n),nc2,wd,ndsr(l,1))
  if ( nad(l,1) == nc2.or.nps /= 1) go to 637
  ne=n+nc2-1
  do 202 i=n,ne
  if ( last == 1)osest=wc(i)
  202 if ( ix(i) == ndsr(l,1)) x(i)=osest
  637 if ( last == 1)go to 2436
  call length(wc(n),x(n),nc2)
  call mstrs(nc2,x(n),wc(n),strss1,strss2)
  strss=strss+strss1
 2436 n=n+nc2
  437 continue
  strss=sqrt(strss/ns)
  return
!
!     row conditional nominal data
!
  509 do 510 l=1,ns
  do 510 i=1,nb
  if ( nsim > 3.and.i <= ncol)go to 503
  call catses(wc(n),x(n),ix(n),nb,wd,ndsr(l,i))
  503 if ( nad(l,i) == nb.or.nps /= 1) go to 610
  if ( nsim > 3.and.i <= ncol)x(n+i-1)=0.0
  ne=n+nb-1
  do 511 j=n,ne
  if ( last == 1)osest=wc(j)
  511 if ( ix(j) == ndsr(l,i)) x(j)=osest
  610 if ( last == 1)go to 1510
  call length(wc(n),x(n),nb)
  call mstrs(nb,x(n),wc(n),strss1,strss2)
  strss=strss+strss1
 1510 n=n+nb
  510 continue
  na=nb
  if ( nsim > 3)na=nb-ncol
  strss=sqrt(strss/(na*ns))
  return
!
!     numerical data
!
  134 rewind ndp
  read(ndp)wd
  rewind ndp
  if ( lin == 1)go to 411
  do 400 i=1,nt
  wdi=wd(i)
  if ( wdi /= big) wd(i)=wdi*wdi
  400 continue
  411 if ( nwc == 1) go to 310
  if ( nwc == 2) go to 512
!
!     unconditional numerical data
!
  call polyf(x,wc,nt,wd,r,alph,ndeg,ncst,nab)
  if ( lin == 1)call lint(x,wc,nt,wd)
  if ( last /= 1)go to 420
  do 415 i=1,nt
  415 if ( x(i) == big)x(i)=wc(i)
  return
  420 call length(wc,x,nt)
  call mstrs(nt,x,wc,strss,strss2)
  strss=sqrt(strss)
  return
!
!     matrix conditional numerical data
!
  310 do 311 l=1,ns
  call polyf(x(n),wc(n),nc2,wd(n),r,alph,ndeg,ncst,nad(l,1))
  if ( lin == 1)call lint(x(n),wc(n),nc2,wd(n))
  if ( last /= 1)go to 300
  do 305 i=1,nc2
  ni=n+i-1
  305 if ( x(ni) == big)x(ni)=wc(ni)
  go to 1311
  300 call length(wc(n),x(n),nc2)
  call mstrs(nc2,x(n),wc(n),strss1,strss2)
  strss=strss+strss1
 1311 n=n+nc2
  311 continue
  strss=sqrt(strss/ns)
  return

!     row conditional numerical data

  512 do 513 l=1,ns
  do 513 i=1,nb
  if ( nsim > 3.and.i <= ncol)go to 514
  if ( nad(l,i) <= 1)go to 514
  call polyf(x(n),wc(n),nb,wd(n),r,alph,ndeg,ncst,nad(l,i))
  if ( lin == 1)call lint(x(n),wc(n),nb,wd(n))
  go to 516
  514 do 515 j=1,nb
  k=j+n-1
  515 x(k)=wc(k)
  516 if ( last /= 1)go to 500
  do 501 j=1,nb
  ni=n+j-1
  501 if ( x(ni) == big)x(ni)=wc(ni)
  go to 1513
  500 call length(wc(n),x(n),nb)
  call mstrs(nb,x(n),wc(n),strss1,strss2)
  strss=strss+strss1
 1513 n=n+nb
  513 continue
  na=nb
  if ( nsim > 3)na=nb-ncol
  strss=sqrt(strss/(na*ns))
  return
end
subroutine driver ( area, eoj )
!
!*****************************************************************************80
!
!! DRIVER controls the flow of the program.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!
!  Discussion:
!
!    The array AREA is a block of memory which is allocated by the
!    main routine, and is subdivided here for use in the subroutines
!
  integer eoj
  dimension area(1)
!
!     integer*2 ix(8000),iy(3000),iz(525),item(55,101)
!     dimension title(20),fmt(20)
!     dimension x(8000),wa(8000),wd(8000)
!     dimension ua(875),u11(36),u12(36),u22(36),r(36),ub1(6),
!    1ub2(6),bk(6),wk(6),wk2(6),wk3(36),wk4(6),xn(6)
!     dimension xx(1600),cfr(280),cfl(280),w(280),tr(40),fk(49)
!     dimension ws(280),ds(1600),zz(225),cv(15),cw(15),zx(280)
!     dimension xeq(160),ndsb(35),ndsr(525),nad(525),fmrr(40)
!     equivalence (ua(1),u11(1)),(ua(37),u12(1)),
!    1(ua(73),u22(1)),(ua(109),r(1)),(ua(145),ub1(1)),(ua(151),ub2(1)),
!    2(ua(157),bk(1)),(ua(163),wk(1)),(ua(169),wk2(1)),(ua(175),wk3(1)),
!    3(ua(211),wk4(1)),(ua(217),xn(1))
!     equivalence (wd(1),item(1,1)),(cfl(1),xeq(1))
!     equivalence (tr(1),fmrr(1)),(zz(1),zx(1))
!----------------------------------------------------------------------

  common /block1/ nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
  common /starts/jx,jwa,jwd,jxx,jds,jcfr,jcfl,jw,jws,jzz,jtr,jfk, &
         jcv,jcw,ju11,ju12,ju22,jub1,jub2,jxn,jphsub,jphsti,jix,jiy,&
         jiz,jndsr,jnad,jftln,jdist,jpijp,jxk
!
  call step1 ( area(jix), area(jx), area(jwa), area(jxx), area(jds), &
    area(jnad), iret )

  if ( iret == 1 ) then
    eoj = 1
    return
  end if

  nparm = nb

  if ( nwe == 1 .or. nwe == 3 ) then
    nparm = nparm + ns
  end if

  if ( nwe > 1 ) then
    nparm = nparm + nb
  end if

  do ijkl = 1, nd

    np = nparm*(ndxp-ijkl)
    if ( np > nt ) go to 7
    if ( nt >= 2.5 * np ) go to 9
    if ( ndtyp < 3 ) go to 9
    write(lout,9910)np,nt
    write(*,9910)np,nt
    call hitr
    go to 9
7   write(lout,9920)np,nt

    write(*,9920)np,nt
    call hitr

    stop

9   continue

    call step2 ( area(jix),area(jiy),area(jiz),area(jx),area(jwa),area(jwa), &
      area(jwd),area(jcfr),area(jcfl),area(jw),area(jtr),area(jfk),area(jxx), &
      area(jws),area(jds),area(jzz),area(jcv),area(jcw),area(jndsr), &
      area(jnad),ijkl )

    call step3 ( area(jix),area(jiy),area(jiz),area(jx),area(jwd), &
      area(jwa),area(jwa),area(jwa),area(ju11),area(ju12), &
      area(ju22),area(jub1),area(jub2),area(jxn),area(jxx), &
      area(jcfl),area(jw),area(jtr),area(jfk),area(jws), &
      area(jndsr),area(jnad),ijkl,ire )

    if ( ire /= 1 ) then

      call step3a ( area(jix), area(jiy), area(jiz), area(jx), &
        area(jwa), area(jwd), area(jxx), area(jcfl), area(jw), &
        area(jws), area(jndsr), area(jnad), area(jphsub), &
        area(jphsti) )

      call step4 ( area(jx), area(jwa), area(jxx), area(jcfl), &
        area(jw), area(jtr), area(jws), &
        area(jdist),area(jds),area(jpijp),area(jxk),area(jtr),ijkl )

    end if

  end do

  return
!
!  warning 9910 has been changed to avoid misunderstanding.
!  the warning is printed, if
!  - the measurement level is ordinal or nominal  .and.
!  - the number of parameters is < 2.5*number of input data.
!
 9910 format(//' alscal warning:',i5,' parameters,',i5,' observations', &
     /8x,'- the number of parameters being computed may be too large', &
     /8x,'- or the number of observations may be too small', &
     /8x,'for reliable results.')

 9920 format(//' alscal fatal error:  computations terminated'/ &
     8x,'the number of parameters (',i5,') exceeds the number', &
     ' of observations(',i5,').')
end
subroutine eigk ( a, value, n, na )
!
!*****************************************************************************80
!
!! EIGK implements Kaiser's JK method for eigenanalysis of a real symmetric matrix.
!
!
!  Discussion:
!
!    This routine finds the eigenvalues and eigenvectors of a real
!    symmetric matrix.  The absolute values of the eigenvalues
!    ordered in descending order of size are returned in value.
!    Normalized eigenvectors are returned in corresponding
!    columns of A.
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     eigen-routine by kaiser (jk method)
!
!  Reference: 
!
!    H F Kaiser,
!    "the eigenvalues of a real symmetric matrix"
!    british computer journal
!    vol 15 no 3 pages 271-273.
!
  integer n
  integer na
!
  real a(na,*)
  real ad
  real d
  double precision halfp
  double precision q
  real value(n)
!
  if ( n == 1 ) then
    value(1) = a(1,1)
    a(1,1) = 1.0
    return
  end if

  if ( n == 2 ) then
    ad=a(1,1)+a(2,2)
    d=sqrt(ad*ad-4.0*(a(1,1)*a(2,2)-a(1,2)*a(2,1)))
    value(1)=(ad+d)*0.5
    value(2)= ad-value(1)
    vl1=1.0/sqrt(a(1,2)**2+(a(1,1)-value(1))**2)
    a(2,1)=(value(1)-a(1,1))*vl1
    vl2=1.0/sqrt(a(1,2)**2+(a(1,1)-value(2))**2)
    a(2,2)=(value(2)-a(1,1))*vl2
    a(1,1)=a(1,2)*vl1
    a(1,2)=a(1,2)*vl2
    return
  end if

  q = sum ( a(1:n,1:n)**2 )
  eps=0.000001*q/n
  nless1=(n-1)
  nn=(n-1)*n/2
  ncount=nn

  116 continue

  do j=1,nless1

    jplus1=j+1

    do k=jplus1,n

      halfp=0.0

      q=0.0
      do i=1,n
        halfp=halfp+a(i,j)*a(i,k)
        q=q+(a(i,j)+a(i,k))*(a(i,j)-a(i,k))
      end do

      absp=dabs(halfp+halfp)

      if ( absp < eps.and.q >= 0.0)go to 106

      absq=dabs(q)

      if ( absp <= absq ) then
        tan=absp/absq
        cos=1.0/sqrt(1.0+tan*tan)
        sin=tan*cos
      else
        ctn=absq/absp
        sin=1.0/sqrt(1.0+ctn*ctn)
        cos=ctn*sin
      end if

      cos=sqrt((1.0+cos)/2.0)
      sin=sin/(cos+cos)

      if ( q < 0.0 ) then
        temp=cos
        cos=sin
        sin=temp
      end if

      if ( halfp < 0.0 ) then
        sin=-sin
      end if

      do i=1,n
        temp=a(i,j)
        a(i,j)=temp*cos+a(i,k)*sin
        a(i,k)=-temp*sin+a(i,k)*cos
      end do

      ncount=nn

      cycle

106   continue

      ncount=ncount-1

      if ( ncount <= 0 )go to 115

    end do

  end do

  go to 116

  115 continue

  do j=1,n

    value(j)=0.0
    do i=1,n
      value(j)=value(j)+a(i,j)*a(i,j)
    end do

    value(j)=sqrt(value(j))

  end do

  do j=1,n
    do i=1,n
      a(i,j)=a(i,j)/value(j)
    end do
  end do

  return
end
subroutine hitr
!
!*****************************************************************************80
!
!! HITR requests that the user hit RETURN, and waits for input.
!
!
!  Discussion:
!
!    For use on "real" computers, these statements have been commented out.
!
!     write ( *, * ) ' '
!     write ( *, * ) 'Hit <return>'
!     read ( *, * ) 

  return
end
subroutine init ( w, cfl, cfr, nb, ns, nd, xx, wa, tr, fk, c, xeq, am, ws, &
  ndx, ndxp )
!
!*****************************************************************************80
!
!! INIT computes the initial configuration estimation in INDSCAL.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     The schonemann-de leeuw method is used.
!
!     Yoshio Takane july 1974
!
  integer nb
  integer ndx
  integer ndxp
  integer ns
!
  real am(ns,*)
  real c(ndx,ndx,*)
  real cfl(nb,*)
  real cfr(nb,*)
  real fk(ndxp,*)
  real w(ns,*)
  real wa(nb,nb,ns)
  real ws(nb,*)
  double precision t1
  double precision t2
  double precision t3
  real tr(*)
  real xeq(ns,*)
  double precision xm
  double precision xs
  real xx(nb,nb)
!
  call cjeig ( xx, cfr, cfl, nb, nd+1, fk, ws, tr, 1, ndxp )
!
!  ws is being passed to cjeig for use as scratch space
!
  do k=1,ns
    do k1=1,nd
      do k2=1,k1
        t1=0.0
        t2=tr(k1)*tr(k2)
        t2=1.0d0/dsqrt(t2)
        do j=1,nb
          t3=0.0
          do i=1,nb
            t3=t3+wa(i,j,k)*cfr(i,k1)
          end do
          t1=t1+t3*cfr(j,k2)
        end do
        t1=t1*t2
        c(k1,k2,k)=t1
        c(k2,k1,k)=t1
      end do
    end do
  end do

  do j=1,nd
    dq=sqrt(tr(j))
    do i=1,nb
      cfr(i,j)=cfr(i,j)*dq
    end do
  end do
!
!  Find feasible weights
!  Forming a matrix

  if ( nd <= 1) go to 301

  do k=1,ns
    t1=0.0
    do j=1,nd
      do i=1,nd
        t1=t1+c(i,j,k)
      end do
    end do
    tr(k)=t1/nb
  end do

  do i=1,ns
    do j=1,i
      t1=0.0
      do k=1,nd
        do l=1,nd
          t1=t1+c(k,l,i)*c(l,k,j)
        end do
      end do
      am(i,j)=t1/nb-tr(i)*tr(j)
      am(j,i)=am(i,j)
    end do
  end do

  call cjeig(am,xeq(1,1),xeq(1,3),ns,2,fk,w,tr,1,ndxp)

  do 18 i=1,nd
  do 18 j=1,i
  t1=0.0
  do 19 k=1,ns
   19 t1=t1+c(i,j,k)*xeq(k,1)
  am(i,j)=t1
   18 am(j,i)=t1
!
!  Find the transformation matrix.
!
  call eigk ( am, tr, nd, ns )
!
!  Find the intitial configuration.
!
  do 20 i=1,nb
    do 20 j=1,nd
      t1=0.0
      do k=1,nd
        t1=t1+cfr(i,k)*am(k,j)
      end do
20 cfl(i,j)=t1
!
!  Find the initial weights
!
  do 22 i=1,nd
  do 22 j=1,i
  t1=0.0
  do 23 k=1,nb
   23 t1=t1+cfl(k,i)*cfl(k,j)
  fk(i,j)=t1
   22 fk(j,i)=t1

  call minv ( fk, nd, ndxp )

  do 21 l=1,ns
  do 24 i=1,nd
  do 24 j=1,i
  t1=0.0
  do 251 k2=1,nb
  t2=0.0
  do 25 k1=1,nb
   25 t2=t2+wa(k1,k2,l)*cfl(k1,i)
  251 t1=t1+t2*cfl(k2,j)
  am(i,j)=t1
   24 am(j,i)=t1
  do 28 i=1,nd
  t1=0.0
  do 26 k1=1,nd
  do 26 k2=1,nd
   26 t1=t1+fk(i,k1)*am(k1,k2)*fk(k2,i)
   28 w(l,i)=t1
   21 continue

!     normalization

  do 57 j=1,nd
  xm=0.0
  xs=0.0
  do 58 i=1,nb
  xm=xm+cfl(i,j)
   58 xs=xs+cfl(i,j)**2
  xm=xm/nb
  xs=xs/nb-xm**2
  do 59 i=1,ns
   59 w(i,j)=w(i,j)*xs
  xs=dsqrt(xs)
  do 60 i=1,nb
   60 cfl(i,j)=(cfl(i,j)-xm)/xs
   57 continue
  return

  301 continue

  do l=1,ns
    w(l,1)=c(1,1,l)
  end do

  return
end
subroutine inner ( cfl, w, x, wb, u11, u12, u22, ub1, ub2, xn, nb, ndim, ns, &
  ndx, nbs, ws )
!
!*****************************************************************************80
!
!! INNER computes stimulus coordinates.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  integer nb
  integer nbs
  integer ndx
  integer ns
!
  real cfl(nb,1)
  double precision t1
  double precision t2
  double precision t3
  double precision t4
  double precision tu11
  double precision tu12
  double precision tu22
  real u11(ndx,1)
  real u12(ndx,1)
  real u22(ndx,1)
  real ub1(1)
  real ub2(1)
  real w(ns,1)
  real wb(nbs,1)
  real ws(nb,1)
  real x(1)
  real xn(1)
!
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /inicon/initx,initw,initws,initxc
!
  nns=nb*ns

  do 15 m=1,nb
    if ( initx == 3.and.m > ncol)go to 15
    if ( initxc == 3.and.m <= ncol)go to 15
    do 13 i=1,ndim
      t1=ws(m,i)
      do 13 j=1,ndim
        t2=t1*ws(m,j)
        tu11=0.0
        tu12=0.0
        tu22=0.0
        do 14 l=1,nb
          if ( l == m)go to 14
          t3=0.0
          do k=1,ns
            t3=t3+w(k,i)*w(k,j)
          end do
          if ( j > i)go to 12
          tu22=tu22+t3*t2
          tu11=tu11+t3*t2*cfl(l,i)*cfl(l,j)
12        tu12=tu12+t3*t2*cfl(l,i)
14      continue
        if ( j > i)go to 13
        u11(i,j)=tu11
        u22(i,j)=tu22
13      u12(i,j)=tu12
    jk=0
    do 39 k=1,ns
      lj=(m-1)*nb
      do 39 j=1,nb
        jk=jk+1
        lj=lj+1
        if ( j == m)go to 153
        t1=wb(lj,k)
        do 40 i=1,ndim
40        t1=t1-w(k,i)*cfl(j,i)**2*ws(m,i)
        x(jk)=t1
        go to 39
153     x(jk)=0.0
39  continue

    do 16 i=1,ndim
      t1=0.0
      t2=0.0
      jk=0
      do 18 k=1,ns
        t3=0.0
        t4=0.0
        do 17 j=1,nb
          jk=jk+1
          t4=t4+x(jk)
17        t3=t3+cfl(j,i)*x(jk)
        t1=t1+t3*w(k,i)
18    t2=t2+t4*w(k,i)
      ub1(i)=t1*ws(m,i)
      ub2(i)=t2*ws(m,i)
16  continue

    if ( nsim < 2) go to 50
    if ( nwe < 2 ) go to 80
    do 84 i=1,ndim
      do 84 j=1,ndim
        tu11=0.0
        tu22=0.0
        tu12=0.0
        do 83 l=1,nb
          if ( l == m)go to 83
          t1=0.0
          do 82 k=1,ns
82          t1=t1+w(k,i)*w(k,j)
          if ( j > i)go to 81
          tu22=tu22+t1*(ws(l,i)*ws(l,j))
          tu11=tu11+t1*(ws(l,i)*ws(l,j))*cfl(l,i)*cfl(l,j)
81        tu12=tu12+t1*(ws(l,i)*ws(l,j))*cfl(l,i)
83      continue
        if ( j > i)go to 84
        u11(i,j)=(u11(i,j)+tu11)*.5d0
        u22(i,j)=(u22(i,j)+tu22)*.5d0
84  u12(i,j)=(u12(i,j)+tu12)*.5d0

80  continue

    jk=nns

    do 59 k=1,ns
      lj=m-nb
      do 59 j=1,nb
        jk=jk+1
        lj=lj+nb
        if ( m == j)go to 253
254     t1=wb(lj,k)
        do 60 i=1,ndim
60        t1=t1-w(k,i)*cfl(j,i)**2*ws(j,i)
        x(jk)=t1
        go to 59
253     x(jk)=0.0
59  continue

    do i=1,ndim
      t1=0.0
      t2=0.0
      jk=nns
      do k=1,ns
        t3=0.0
        t4=0.0
        do j=1,nb
          jk=jk+1
          t3=t3+x(jk)*ws(j,i)*cfl(j,i)
          t4=t4+x(jk)*ws(j,i)
        end do
        t1=t1+w(k,i)*t3
        t2=t2+w(k,i)*t4
      end do
      ub1(i)=(ub1(i)+t1)*.5d0
      ub2(i)=(ub2(i)+t2)*.5d0
    end do

50  continue

    do 68 i=1,ndim
      do 86 j=1,ndim
        if ( j > i)go to 86
        u11(i,j)=4.0*u11(i,j)
        u11(j,i)=u11(i,j)
        u22(j,i)=u22(i,j)
86    u12(i,j)=-2.*u12(i,j)
68  ub1(i)=-2.0*ub1(i)

    call coef(u11,u12,u22,ub1,ub2,ndim,cfl,m,xn,nb,ndx)

15 continue

  return
end
subroutine inswm ( ds, cfl, cfr, ws, xx, tr, cv, cw, fk, zz, c, nadct, nb, &
  ndim, nd, ndx, ndxp, ns )
!
!*****************************************************************************80
!
!! INSWM computes the initial configuration for the stimulus weighted model.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  integer nb
  integer nd
  integer ndx
  integer ndxp
!
  real c(ndx,ndx,1)
  real cfl(nb,1)
  real cfr(nb,1)
  real fk(ndxp,1)
  double precision t
  double precision tt
  real ws(nb,1)
  real xx(nb,nb)
  real zz(nd,1)

  dimension ds(nb,nb),tr(1),cv(1),cw(1)
!
  if ( nadct /= 1 ) then

    t=0.0
    do j=1,nb
      do i=1,nb
        t=t+ds(i,j)**2
      end do
    end do

    t=t/((nb*ns)**2)

    do j=1,nb
      do i=1,nb
        ds(i,j)=ds(i,j)/t
      end do
    end do

    nadct=1

  end if

  do j=1,ndim
    cfr(1,j)=cw(j)
  end do

  do 20 i=1,nb
  do 21 j=1,nb
  do 22 k=1,ndim
   22 xx(j,k)=(cfl(i,k)-cfl(j,k))**2*cfr(1,k)
  if ( ndim == 1) go to 21
  n=ndim
  do 23 k=2,ndim
  nm2=k-1
  do 23 kk=1,nm2
  n=n+1
   23 xx(j,n)=(cfl(i,kk)-cfl(j,kk))*(cfl(i,k)-cfl(j,k))* &
    2.0*sqrt(cfr(1,k)*cfr(1,kk))
   21 continue
  do 24 k=1,nd
  do 24 j=1,nd
  zz(j,k)=0.0
  do 24 kk=1,nb
   24 zz(j,k)=zz(j,k)+xx(kk,j)*xx(kk,k)

  call minv(zz,nd,nd)

  do 25 j=1,nd
  cw(j)=0.0
  do 25 k=1,nb
   25 cw(j)=cw(j)+xx(k,j)*ds(i,k)
  do 26 j=1,nd
  cv(j)=0.0
  do 26 k=1,nd
   26 cv(j)=cv(j)+zz(j,k)*cw(k)
  do 27 j=1,ndim
   27 c(j,j,i)=cv(j)
  if ( ndim == 1) go to 151
  n=ndim
  do 28 k=2,ndim
  nm2=k-1
  do 28 kk=1,nm2
  n=n+1
  c(k,kk,i)=cv(n)
   28 c(kk,k,i)=cv(n)
  151 continue
   20 continue
  do 29 k=1,nb
  tr(k)=0.0
  do 30 j=1,ndim
  do 30 i=1,ndim
   30 tr(k)=tr(k)+c(i,j,k)
   29 tr(k)=tr(k)/nb
  do 31 j=1,nb
  do 31 i=1,nb
  t=0.0
  do 32 k=1,ndim
  do 32 l=1,ndim
   32 t=t+c(k,l,i)*c(l,k,j)
   31 xx(i,j)=t/nb-tr(i)*tr(j)

  if ( ndim == 1) go to 130
  call cjeig(xx,ws,cfr,nb,2,fk,zz,tr,1,ndxp)
  go to 131
  130 do 132 k=1,nb
  132 ws(k,1)=tr(k)
  131 continue

  do 33 j=1,ndim
  do 33 i=1,ndim
  xx(i,j)=0.0
  do 133 k=1,nb
  133 xx(i,j)=xx(i,j)+c(i,j,k)*ws(k,1)
   33 continue

  call eigk ( xx, tr, ndim, nb )

  do 60 j=1,ndim
  do 64 i=1,nb
  t=0.0
  do 63 k=1,ndim
   63 t=t+cfl(i,k)*xx(k,j)
   64 cfr(i,j)=t
   60 continue
  do 35 j=1,ndim
  do 35 i=1,nb
  cfl(i,j)=cfr(i,j)
  t=0.0
  do kk=1,ndim
    tt=0.0
    do k=1,ndim
      tt=tt+xx(k,j)*c(k,kk,i)
    end do
    t=t+tt*xx(kk,j)
  end do
   35 ws(i,j)=t

  return
end
subroutine length ( dist, disp, n )
!
!*****************************************************************************80
!
!! LENGTH normalizes the length of disparities.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     normalize length of disparities so that they optimize
!        normalized sstress insted of raw stress.
!        adjustment is based only on distances which correspond to
!        active data.  adjustment is only made to disparities which
!        correspond to active data.  type of adjustment depends on
!        stress formula being optimized.
!

  dimension dist(1),disp(1)
  double precision ssq,spd,ratio,distv,prod,distm,dispm

  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol

  if ( nsim > 3)go to 3

!     adjustment for stress formula 1

  ssq=0.0
  spd=0.0
  do 1 i=1,n
  if ( disp(i) == big)go to 1
  ssq=ssq+dist(i)**2
  spd=spd+dist(i)*disp(i)
    1 continue
  if ( spd == 0.0)return
  ratio=ssq/spd
  do 2 i=1,n
  if ( disp(i) == big)go to 2
  disp(i)=disp(i)*ratio
    2 continue
  return

!     adjustment for stress formula 2

    3 distm=0.0
  dispm=0.0
  nact=0
  do 4 i=1,n
  if ( disp(i) == big)go to 4
  distm=distm+dist(i)
  dispm=dispm+disp(i)
  nact=nact+1
    4 continue
  distm=distm/nact
  dispm=dispm/nact
  prod=0.0
  distv=0.0
  do 5 i=1,n
  if ( disp(i) == big)go to 5
  distv=distv+(dist(i)-distm)**2
  prod=prod+(dist(i)-distm)*(disp(i)-dispm)
    5 continue
  if ( prod == 0.0)return
  ratio=distv/prod
  do 6 i=1,n
  if ( disp(i) == big)go to 6
  disp(i)=disp(i)*ratio+distm*(1.0-ratio)
    6 continue

  return
end
subroutine lint ( disp, dist, n, obs )
!
!*****************************************************************************80
!
!! LINT performs constrained least squares linear regression.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     routine to perform constrained least squares linear regression.
!        it finds the best linear transformation for quantitative
!        (interval level) data under the constraint that all
!        disparities (predictions) be non-negative.
!        the newton-raphson method is used.

!   original version by Yoshio Takane and Forrest Young
!   rewritten by Rostyslaw Lewyckyj   march 1977
!   modified according to suggestions by Yoshio Takane june 1977
!   rewritten in double precision by Yoshio Takane september 1978
!   errors corrected by Yoshio Takane and Forrest Young, january 1981
!   adapted for stand alone version of alscal, july 1982
!
  implicit double precision (a-h,o-z)
  real disp(1),dist(1),obs(1)
  real prtmsg
!
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata

  double precision cut,stmin
  integer debug,icnstr,noulb
  common /prmblk/cut,stmin,debug,icnstr,noulb
  data prtmsg/0.0/
  data big/9.0e20/
!---the following statement is sas only
!     iopt19=iopt(19)
  fmin=obs(1)
  do 10 k=1,n
   10 if ( fmin > obs(k))fmin=obs(k)
  dh=0.0
  odh=0.0
  o1=0.0
  o2=0.0
  o3=0.0
  o4=0.0
  dc=0.0
  d0=0.0
  d1=0.0
  d2=0.0

! compute initial estimates by linear regression

  do 22 k=1,n
  if ( obs(k) == big) go to 22
  di=sqrt(dist(k))
  dh=dh+di
  ob1=obs(k)
  odh=odh+ob1*di
  di=dist(k)
  o1=o1+ob1
  dc=dc+1.0
  d0=d0+di
  d1=d1+ob1*di
  ob=ob1*ob1
  o2=o2+ob
  d2=d2+ob*di
  ob=ob*ob1
  o3=o3+ob
  o4=o4+ob*ob1
   22 continue
  ao=(odh*dc-dh*o1)/(o2*dc-o1*o1)
  bo=(dh-ao*o1)/dc
  fo=ao*ao*ao*ao*o4+6.0*ao*ao*bo*bo*o2+dc*bo*bo*bo*bo-2.0*ao*ao*d2 &
    -4.0*ao*bo*d1-2.0*bo*bo*d0+4.0*ao*ao*ao*bo*o3+4.0*ao*bo*bo*bo*o1
!
! perform unconstrained newton-raphson iterations
!
    4 ll=0
    1 ll=ll+1
  if ( ll > 50)go to 2
  g1=((o4*ao+3.0*bo*o3)*ao+3.0*bo*bo*o2-d2)*ao+bo*(bo*bo*o1-d1)
  g2=((dc*bo+3.0*ao*o1)*bo+3.0*ao*ao*o2-d0)*bo+ao*(ao*ao*o3-d1)
  h11=d2-3.0*ao*ao*o4-6.0*ao*bo*o3-3.0*bo*bo*o2
  h22=d0-3.0*ao*ao*o2-6.0*ao*bo*o1-3.0*bo*bo*dc
  h12=d1-3.0*ao*ao*o3-6.0*ao*bo*o2-3.0*bo*bo*o1
  det=h11*h22-h12*h12
  det=1.0/det
  h12=-h12*det
  p=h11
  h11=h22*det
  h22=p*det
  s1=h11*g1+h12*g2
  s2=h12*g1+h22*g2
  step=1.0
  lll=0
    3 lll=lll+1
  if ( lll > 20) go to 2
  an=ao+step*s1
  bn=bo+step*s2
  fn=an*an*an*an*o4+6.0*an*an*bn*bn*o2+dc*bn*bn*bn*bn-2.0*an*an*d2 & 
    -4.0*an*bn*d1-2.0*bn*bn*d0+4.0*an*an*an*bn*o3+4.0*an*bn*bn*bn*o1
  if ( fn < fo) go to 40
  if ( fn == fo) go to 2
  step=step*0.5
  go to 3
   40 gn=h11*g1*g1+2.0*h12*g1*g2+h22*g2*g2
  gn=dsqrt(-gn)
  fo=fn
  if ( gn <= 0.0001) go to 2
  ao=an
  bo=bn
  go to 1

! test for negative slope.

    2 if ( an < 0.0) go to 50

! test for negative estimates

  222 if ( an*fmin+bn < 0.0) go to 30

! come here when slope and estimates are positive

  do 12 k=1,n
  if ( obs(k) == big)go to 13
  disp(k)=(an*obs(k)+bn)*(an*obs(k)+bn)
  go to 12
   13 disp(k)=big
   12 continue
!---the following statement is replaced by the next
!     if ( iopt19 == 1)write(lout,4738)an,bn
  if ( debug > 0)write(lout,4738)an,bn
 4738 format(' lint: positive est+++slope=',e15.5,'intercept=',e15.5)
  return

! come here when the slope is negative.  set the intercept to zero
! (ratio level), print warning, and solve for transformation.

   50 fmin=0.0
  if ( prtmsg == 1)go to 30
  prtmsg=1
  write(lout,9998)
  write(lout,9995)

! come here when negative estimates do exist and impose constraint

   30 t1=0.0
  t2=0.0
  do 36 k=1,n
  if ( obs(k) == big)go to 36
  s=(obs(k)-fmin)**2
  t1=t1+s*dist(k)
  t2=t2+s**2
   36 continue
  a=t1/t2
  do 37 k=1,n
  if ( obs(k) == big)go to 38
  disp(k)=a*(obs(k)-fmin)**2
  go to 37
   38 disp(k)=big
   37 continue
!---the following statement replaced by the next 7/14/82
!     if ( iopt19 == 1)write(lout,4748)a,fmin
  if ( debug > 0)write(lout,4748)a,fmin
 4748 format(' lint: negative est---slope=',e15.5,'intercept=',e15.5)
  return
 9995 format(t18,'the measurement level has been changed to ratio on'/ &
   ,t18,'this iteration and will be changed back to interval on'/ &
   ,t18,'the next.  this may result in divergence.  if so,'/ &
   ,t18,'reanalyze your data with a different measurement level.'/ &
   ,t18,'in all cases you should check the options of the'/ &
   ,t18,'analysis, particularly the similar option.')

 9998 format(/' alscal warning: a linear transformation of your data has' // &
    'negative slope.')
end
subroutine minv ( a, n, na )
!
!*****************************************************************************80
!
!! MINV computes the inverse of a symmetric matrix.
!
  dimension a(1)
  dimension l(30)
  dimension m(30)
!
  if ( n  ==  1 ) then

    a(1) = 1.0 / a(1)
    return

  else if ( n  ==  2 ) then

    c=1.0/(a(1)*a(na+2)-a(na+1)*a(2))
    w=a(1)
    a(1)=a(na+2)*c
    a(na+2)=w*c
    a(na+1)=-a(na+1)*c
    a(2)=-a(2)*c
    return

  end if

!  pack the matrix

  k=n
  i1=na+1
  i2=n*na

  do i3=i1,i2,na
    i4=i3+n-1
    do i=i3,i4
      k=k+1
      a(k)=a(i)
    end do
  end do

!  search for largest element

  nk=-n
  do 80 k=1,n
  nk=nk+n
  l(k)=k
  m(k)=k
  kk=nk+k
  biga=a(kk)
  do 20 j=k,n
  iz=n*(j-1)
  do 20 i=k,n
  ij=iz+i
  if ( abs(biga) < abs(a(ij)))go to 20
  biga=a(ij)
  l(k)=i
  m(k)=j
   20 continue

!  interchange rows

  j=l(k)
  if ( j <= k)go to 35
  ki=k-n
  do 30 i=1,n
  ki=ki+n
  hold=-a(ki)
  ji=ki-k+j
  a(ki)=a(ji)
   30 a(ji) =hold

!  interchange columns

   35 i=m(k)
  if ( i <= k)go to 45
  jp=n*(i-1)
  do 40 j=1,n
  jk=nk+j
  ji=jp+j
  hold=-a(jk)
  a(jk)=a(ji)
   40 a(ji) =hold

!  divide column by minus pivot (value of pivot element is
!  contained in biga)

   45 if ( abs(biga) > 1.0e-20)go to 48
  900 print 901
  901 format(/' alscal fatal error:  inverse of a singular matrix attempted.')
!
!  call errtra
!
  stop
   48 do 55 i=1,n
  if ( i == k)go to 55
  ik=nk+i
  a(ik)=a(ik)/(-biga)
   55 continue

!  reduce matrix

  do 65 i=1,n
  ik=nk+i
  hold=a(ik)
  ij=i-n
  do 65 j=1,n
  ij=ij+n
  if ( i == k)go to 65
  if ( j == k)go to 65
  kj=ij-i+k
  a(ij)=hold*a(kj)+a(ij)
   65 continue

!  divide row by pivot

  kj=k-n
  do 75 j=1,n
  kj=kj+n
  if ( j /= k)a(kj)=a(kj)/biga
   75 continue

!  product of pivots
!  replace pivot by reciprocal

  a(kk)=1.0/biga
   80 continue

!  final row and column interchange

  k=n
  100 k=(k-1)
  if ( k <= 0)go to 150
  105 i=l(k)
  if ( i <= k)go to 120
  108 jq=n*(k-1)
  jr=n*(i-1)
  do 110 j=1,n
  jk=jq+j
  hold=a(jk)
  ji=jr+j
  a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
  if ( j <= k)go to 100
  125 ki=k-n
  do 130 i=1,n
  ki=ki+n
  hold=a(ki)
  ji=ki-k+j
  a(ki)=-a(ji)
  130 a(ji) =hold
  go to 100

!  unpack the matrix

  150 i3=n*(na+1)+1
  i4=n*(n+1)+1

  do j=2,n
    i3=i3-na
    i4=i4-n
    do i=1,n
      a(i3-i)=a(i4-i)
    end do
  end do

  return
end
subroutine mstrs ( n, disp, dist, strss1, resid )
!
!*****************************************************************************80
!
!! MSTRS calculates the stress of two vectors and estimates missing data.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     stress 1 or 2 is calculated (normalized by disparities, which is
!        the same as by distances).  the residual variance (1-rsq)
!        is calculated for normalization of derived weights.
!
!     estimates of missing data are their distances
!

  dimension disp(1),dist(1)
  double precision sv4,sv4sq,ssq,sxy,s4,s4sq
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strssb,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
!
  nact=0
  sxy=0.0
  s4=0.0
  s4sq=0.0
  sv4=0.0
  sv4sq=0.0
  ssq = 0

  do i = 1, n

    if ( disp(i) /= big ) then
      s4=s4+dist(i)
      sv4=sv4+disp(i)
      sv4sq=sv4sq+disp(i)*disp(i)
      s4sq=s4sq+dist(i)*dist(i)
      sxy=sxy+dist(i)*disp(i)
      ssq=ssq+(disp(i)-dist(i))**2
      nact=nact+1
    else
      disp(i)=dist(i)
    end if

  end do

  if ( nact <= 1) then
    strss1=0.0
    resid=0.0
    return
  end if

  if ( nsim < 4 ) go to 118

  dispm=sv4/nact

  svar=0.0
  do i=1,n
    if ( disp(i) /= big ) then
      svar=svar+(disp(i)-dispm)**2
    end if
  end do

  strss1=ssq/svar
  go to 119
  118 strss1=ssq/sv4sq
  119 resid=1.0-(sxy-s4*sv4/nact)**2/((s4sq-s4*s4/nact)*(sv4sq-sv4*sv4/nact))
  return
end
subroutine normw ( w, n, nd, phi, phirow, iconfl )
!
!*****************************************************************************80
!
!! NORMW normalizes weight matrices for output.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  dimension w(n,1),phirow(1)
  double precision sumsq

!     normalize either stimulus or subject weights so that their
!     length corresponds to the proportion  of variance in the
!     optimally scaled data accounted for by the model.

  sumsq=0.0

  do 3 i=1,n

    do j=1,nd
      sumsq=sumsq+w(i,j)**2
    end do

    if ( iconfl == 1) go to 3
!
!     come here when conditionality is such that each row of weights
!     is normalized separately.
!
    sumsq=dsqrt(phirow(i)/sumsq)
    do j=1,nd
      w(i,j)=sumsq*w(i,j)
    end do

    sumsq=0.0

    3 continue

  if ( iconfl == 0) return

!     come here when conditionality is such that rows are
!     normalized jointly.

  sumsq=dsqrt(n*phi/sumsq)

  do j=1,nd
    do i=1,n
      w(i,j)=sumsq*w(i,j)
    end do
  end do

  return
end
subroutine normx ( x, w, nb, ns, nd, nwe )
!
!*****************************************************************************80
!
!! NORMX normalizes the configuration for printing at end of job.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  dimension x(nb,1),w(ns,1)
  double precision sumsq,sum

  sumsq=0.0
  do 4 j=1,nd

!     center configuration at centroid

  sum=0.0
  do i=1,nb
    sum=sum+x(i,j)
  end do

  sum=sum/nb
  do i=1,nb
    x(i,j)=x(i,j)-sum
    sumsq=sumsq+x(i,j)**2
  end do

  if (nwe == 0) go to 4
!
!     for all of the weighted models normalize so that the length of
!     each stimulus dimension is unity
!
  sumsq=dsqrt(nb/sumsq)
  do i=1,nb
    x(i,j)=sumsq*x(i,j)
  end do

  sumsq=1.0/(sumsq**2)

  do i=1,ns
    w(i,j)=sumsq*w(i,j)
  end do

  sumsq=0.0
    4 continue
  if ( nwe > 0) return
!
!  For the unweighted model normalize length of the dimensions so
!  that their average length is unity.
!
  sumsq=dsqrt(nb*nd/sumsq)

  x(1:nb,1:nd) = sumsq * x(1:nb,1:nd)

  return
end
subroutine outa ( x, nb, nfl )
!
!*****************************************************************************80
!
!! OUTA prints out asymmetric and rectangular matrices.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  dimension  x(nb,nb,1)
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata

  ns=abs(nfl)
  na=1
  if ( nsim > 3)na=ncol+1
  nc=nb
  if ( nsim > 3)nc=ncol

  do m=1,ns

    if ( nfl > 0)write(lout,100)m
  100   format(//' matrix',i4)

    do i1=1,nc,10

      i2 = min ( nc, i1+9 )
      write(lout,200)(k,k=i1,i2)
  200     format(/ 13x,10i10)

      n=0
      do l=na,nb
        n=n+1
        if ( nfl < 0)then
          write(lout,201)n,(x(l,k,m),k=i1,i2)
        else
          write(lout,201)n,(x(k,l,m),k=i1,i2)
        end if
      end do

    end do

  end do

  201 format(i16,10f10.3)
  return
end
subroutine outs ( x, nb, ns, nc, xx )
!
!*****************************************************************************80
!
!! OUTS prints out symmetric matrices.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  integer n
  real x(1)
  real xx(nb,1)
!
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
!
  n = 0
  do i = 1, nb
    xx(i,i) = 0.0
  end do

  do m = 1, ns

    if ( nc  >  0 ) then

      write(lout,100) m
  100     format(//'   matrix',i5)

      do i = 2, nb
        i1=i-1
        do j=1,i1
          n=n+1
          xx(i,j)=x(n)
          xx(j,i)=x(n)
        end do
      end do

    end if

    do i1 = 1, nb, 10
      i2 = min(nb,i1+9)
      write(lout,200)(k,k=i1,i2)
      do l = i1, nb
        i3=min ( l, i2 )
        write(lout,201)l,(xx(l,k),k=i1,i3)
      end do
    end do

  200   format(/ 13x,10i10)
  201   format(i16,10f10.3)

  end do

  return
end
subroutine page ( k )
!
!*******************************************************************************
!
!! PAGE writes a form feed.
!
  write ( k, * ) ' '
  write ( k, * ) ' '
  write ( k, * ) ' '
  write ( k, * ) ' '

!     write ( k, '(a1)' ) char(12)

  return
end
subroutine pddisp ( o, n, isim, ifl, sumsq )
!
!*******************************************************************************
!
!! PDDISP prepares disparities for PDSCAL.
!
!
!  Parameters:
!
!      o:    disparities (n by n)
!      isim: true=similarity
!      ifl:  true=asymmetric
!      sumsq:sums of squares of disparities after normalization
!
  integer n
!
  logical ifl
  logical isim
  integer na
  integer nb
  real o(n,n)
!
   nb = n
   na = 1
   if ( .not. ifl ) then
     na=2
   end if

19     if ( .not.isim)go to 29

! if similarity make dissimilarity

   xmax = 0
   do i=na,n
     if ( .not.ifl)nb=i-1
     do j=1,nb
       if ( xmax <= o(i,j) ) then
         xmax=o(i,j)
       end if
     end do
   end do

   do i=na,n
     if ( .not.ifl)nb=i-1
     do j=1,nb
       o(i,j)=xmax-o(i,j)+1
     end do
   end do
!
! square all entries

 29    do 30 i=na,n
   if ( .not.ifl)nb=i-1
   do 30 j=1,nb
   if ( i == j)go to 30
   o(i,j)=o(i,j)**2
 30    continue

   if ( 1 == 1 ) go to 61

!      normalize disparities so that the sum of
!      squared disparities equals number of elements

  total=0.0
  do i=na,n
    if ( .not.ifl)nb=i-1
    do j=1,nb
      if ( i /= j ) then
        total=total+o(i,j)**2
      end if
    end do
  end do

  if ( ifl)xnum=n*(n-1)
  if ( .not.ifl)xnum=n*(n-1)/2
  xnum=sqrt(xnum)
  sqtot=sqrt(total)
  do i=na,n
    if ( .not.ifl)nb=i-1
    do j=1,nb
      if ( i /= j ) then
        o(i,j)=o(i,j)*xnum/sqtot
      end if
    end do
  end do

!       fill in upper half of symmetric data
!       (only needed for initialization routine)

61    continue
 
  if ( .not. ifl ) then

    do i=1,n
      o(i,i)=0.0
      do j=i,n
        o(i,j)=o(j,i)
      end do
    end do
  end if
!
!  calculate s-stress denominator
!
  sumsq=0.0

  do i=na,n
    if ( .not.ifl)nb=i-1
    do j=1,nb
      sumsq=sumsq+o(i,j)**2
    end do
  end do

  return
end
subroutine pddist ( x, a, n, r, t, d )
!
!*******************************************************************************
!
!! PDDIST calculates principal direction distances.
!
!
!  Author:
!
!    Forrest W. Young and Cynthia Null
!
!  Parameters:
!
!    x:multivariate data(common space)
!    a:principal direction weights
!    n:number of observations(stimuli)
!    r:number of variables(attributes)
!    t:number of principal directions
!    d:distances in the principal directions space
!
  integer n
!
  real a(6,6)
  real d(n,n)
  integer r
  integer t
  real x(n,6)
!
  do i=2,n
    do j=1,i-1
      d(i,j)=0.0
      do ia=1,r
        pija=x(i,ia)-x(j,ia)
        do ib=1,r
          prod=(x(i,ib)-x(j,ib))*pija
          do ic=1,t
            d(i,j)=d(i,j)+prod*a(ia,ic)*a(ib,ic)
          end do
        end do
      end do
    end do
  end do

  return
end
subroutine pdinwt ( x, o, xk, a, b, n, r, t, row, xtx, xtxk, ifl )
!
!*******************************************************************************
!
!! PDINWT obtains an initial estimate of the principal direction weights.
!
!
  real a(6,6)
  integer r
  integer t

  real xtx(6,6),xtxk(6,6)
  real o(n,n),x(n,6),xk(n,6),row(n),b(n,n)
  logical ifl
!
!  compute scalar products from this subjects dissimilarities
!
  do i = 1, n
    row(i) = sum ( o(i,1:n) ) / real ( n )
  end do

  grand = sum ( row(1:n) ) / real ( n )

  do i = 1, n
    do j = 1, n
      b(i,j)=-0.5*(o(i,j)-row(i)-row(j)+grand)
    end do
  end do
!
!  symmetrize asymmetric scalar products
!
  if ( ifl ) then

    do i=2,n
      do j=1,i-1
        b(i,j)=(b(i,j)+b(j,i))*.5
        b(j,i)=b(i,j)
      end do
    end do

  end if
!
!  Compute unconstrained personal space coordinates
!  by the classical torgerson procedure.

  call eigk ( b, row, n, n )

  do j=1,t
    row(j)=sqrt(row(j))
    do i=1,n
      xk(i,j)=b(i,j)*row(j)
    end do
  end do
!
! compute linear least squares approximation to
! principal direction weights by regression
!
  do i=1,r
    do j=1,r
      xtx(i,j)=0.0
      do k=1,n
        xtx(i,j)=xtx(i,j)+x(k,i)*x(k,j)
      end do
    end do
  end do

  call minv(xtx,r,6)

  do i=1,r
    do j=1,t
      xtxk(i,j)=0.0
      do k=1,n
        xtxk(i,j)=xtxk(i,j)+x(k,i)*xk(k,j)
      end do
    end do
  end do

  do i=1,r
    do j=1,t
      a(i,j)=0.0
      do k=1,r
        a(i,j)=a(i,j)+xtx(i,k)*xtxk(k,j)
      end do
    end do
  end do

  return
end
subroutine pdmain ( x, o, n, r, t, isub, isym, iplot, iout, d, uk, pijp, xk, &
  row )
!
!*******************************************************************************
!
!! PDMAIN is the main subroutine for principal direction scaling.
!
!
!                             subroutines for
!                      principal  directions scaling
!
!         pdscal82         evolved from pdscal80          dec 1981
!                                 14dec81
!
!                           final change 23may83
!
!               an algorithm to find a subspace  of a space
!               such that the squared distances between the
!               stimuli in the subspace are a least squares
!               fit to a  matrix  of  dissimilarities  data
!
!         copyright 1980  by Forrest W. Young and  cynthia h. null
!
!         the details of the analysis and algorithm  are given in:
!         Young, f.w. principal directions scaling  (notes 1,2 and
!         3) psychometric laboratory, university of north carolina
!         1979 (xeroxed rough drafts). chapel hill, n.c. 27514 usa
!
! input structure
! ---------------
!
! x     =  common stimulus space (n by r)
! o     =  disparities - scaled data (n by n) must be square
! n     =  number of stimuli
! r     =  number of dimensions
! t     =  number of principal directions (t <= r)
! isym  <2 symmetric disparities
!       >1 asymmetric disparities
! iplot =0 for no plots
!       =1 to plot the subject space
! iout  =  printer logical unit number
!
!          the remaining input variables are all internal to
!          this routine and are in the subroutine statement
!          so that space can be dynamically allocated for them
!
! fixed values of internal variables
! ----------------------------------
!
! itmax =50    maximum number of iterations
! crit  =.001  upper limit on coefficient change to end iterations
! ifull =0     do not print debugging iteration history
!       =1     do print debugging iteration history

  parameter ( itmax = 50 )
  parameter ( icrit = 0.0005 )

  real ck(6,6),a(6,6),xtx(6,6),xtxk(6,6),wk(6,6)
  real o(n,n),d(n,n),pijp(n,n)
  real x(n,6),xk(n,6),uk(n,6),row(n)
  equivalence (a(1,1),ck(1,1))
  logical*1 xceed(6,6)
  integer r,t
  logical ifl,isim
  double precision cut,stmin
  integer debug,icnstr,noulb
!
  common /prmblk/cut,stmin,debug,icnstr,noulb
!
  ifl=.false.
  if ( isym > 1)ifl=.true.
  isim=.false.
  if ( isym == 1.or.isym == 3.or.isym == 5)isim=.true.
  ifull=debug
!
!  Massage and print disparities
!
  call pddisp(o,n,isim,ifl,sumsq)
!
!  Compute initial subspace coefficients
!
  call pdinwt(x,o,xk,a,d,n,r,t,row,xtx,xtxk,ifl)
!
!  Compute optimal subspace coefficients
!
  call pdscal(a,x,o,d,n,r,t,ifl,iout,itmax,crit,pijp,sumsq,xceed,ifull)
!
!  Calculate, and print Young's principal directions solution,
!  including the weight matrix (wk),
!  principal direction coefficients (ck),
!  and the normalized orthogonal personal subspace (xk).
!
  call pdyung(ck,wk,xk,x,n,r,t,row)
  call pdouty(ck,wk,xk,n,r,t,isub,iout)
!
!  Calculate, and print tucker's oblique axes solution
!  including the correlations between the oblique axes (wk),
!  weights aplied to each oblique axis (diag(wk)),
!  and the projections onto the oblique axes (uk).
!
  call pdtuck(wk,x,uk,n,r)
  call pdoutt(wk,uk,n,r,isub,iout)
!
!  Plot Young's principal directions space
!
  if ( iplot /= 0) then
    call pdplot(xk,n,t,iout,n)
  end if

  return
end
subroutine pdoutt ( wk, uk, n, r, isub, iout )
!
!*******************************************************************************
!
!! PDOUTT prints Tucker's oblique decomposition of a subject's matrix of principal weights.
!
!
!  Parameters:
!
!     wk   axis correlations (off diagonal) and weights (on diagonal)
!     uk   axis projections
!     n    stimuli
!     r    variables/axes
!     isub subject number
!     iout output unit number
!
  character ( len = 80 ) fmt
  integer r
  character ( len = 80 ) title
  real uk(n,6)
  real wk(6,6)
!
!-----the following statement added 8/13/82
  common /block3/title,fmt
  call page ( iout )
  write(iout,101)title,isub,(i,i=1,r)
101   format(a80//' general Euclidean model:  subject',i4/ &
     ' tucker''s oblique decomposition of the generalized weight matrix.'/ &
     /' correlations between oblique axes'/'   axis',14x,'axis'/3x,(7i10))
  i=1
  one=1.0
  write(iout,100)i,one
100   format(i6,7f10.4)

  do i = 2, r
    write(iout,100)i,(wk(i,j),j=1,i-1),one
  end do

  write(iout,102)isub,(i,i=1,r)
102   format(/' weights on oblique axes for subject',i3/ &
    ' subject',13x,'axis'/3x,(7i10))
  write(iout,100)isub,(wk(i,i),i=1,r)
  write(iout,103)isub,(i,i=1,r)
103   format(/' oblique space for subject',i3/' stimulus',12x,'axis'/3x,(7i10))

  do i=1,n
    write(iout,100)i,(uk(i,j),j=1,r)
  end do

  return
end
subroutine pdouty ( ck, wk, xk, n, r, t, isub, iout )
!
!*******************************************************************************
!
!! PDOUTY prints out CK, WK and XK.
!
!
!  Parameters:
!
!       ck      principal direction coefficients
!       wk      principal direction weights (transformation)
!       xk      principal direction subspace
!       n       number of stimuli
!       r       number of principal directions
!       t       number of variables
!       isub    subject number
!       iout    output unit number
!
  character ( len = 80 ) fmt
  integer r
  integer t
  character ( len = 80 ) title
  real xk(n,6),ck(6,6),wk(6,6)
!
  common /block3/title,fmt
  write(iout,820)isub,(i,i=1,r)
820   format(/' generalized weight matrix for subject',i3/ &
    ' dimension',8x,'dimension'/3x,(7i10))

  do i=1,r
    write(iout,823)i,(wk(i,j),j=1,r)
  end do

  call page ( iout )
  write(iout,8820)title,isub
8820  format(a80//' general Euclidean model:  subject',i4/ &
    1x,' Young''s orthogonal principal directions', &
     ' decomposition of the generalized weight matrix.')
100   write(iout,822)(i,i=1,t)
822   format(/' orthogonal principal direction coefficients'/ &
    ' dimension        direction'/3x,(7i10))

  do j=1,r
    write(iout,823)j,(ck(j,i),i=1,t)
823     format(i6,(7f10.4))
  end do

  write(iout,824)isub,(i,i=1,t)
824   format(/' personal subspace for subject',i3/ &
     ' stimulus         direction'/3x,(7i10))

  do i=1,n
    write(iout,823)i,(xk(i,j),j=1,t)
  end do

  return
end
subroutine pdplot ( xk, n, t, iout, nmax )
!
!*******************************************************************************
!
!! PDPLOT plots each pair of directions of XK.
!
!
!  Parameters:
!
!      xk:principal direction space
!      n:number of stimuli
!      t:number of principal directions
!      nmax:number number of stimuli allowed
!      iout:out unit number
!
!      plot xk
!
  character ( len = 80 ) fmt
  integer t
  character ( len = 80 ) title
  real xk(n,6)
!
  common /block3/title,fmt
  common /page_common/ nlines
!
  if ( t /= 1 ) then

    do j = 2, t

      do i = 1, j-1

        call page ( iout )
        write(iout,798)title,i,j
798         format(a80//' plot of principal direction',i2, &
            ' (horizontal) vs.',1i2,' (vertical).')
        call plotr(xk(1,i),xk(1,j),2.5,2.5,-2.5,-2.5,n,iout,2,nmax)
      end do
    end do

  else

    call page ( iout )
    write(iout,799)title
799     format(a80//' plot of first principal direction')
    call plotr(xk(1,1),xk(1,1),2.5,2.5,-2.5,-2.5,n,iout,2,nmax)

  end if

  return
end
subroutine pdscal ( a, x, o, d, n, r, t, ifl, iout, itmax, crit, pijp, sumsq, &
  xceed, ifull )
!
!*******************************************************************************
!
!! PDSCAL computes principal direction weights.
!
!
!  Parameters:
!
!         x:multivariate data(common space)
!         o:dissimilarity data(or disparities)
!         d:squared Euclidean distances in the principal
!           directions space
!         n:number of observations (stimuli)
!         r:number of variables (attributes)
!         t:number of principal directions
!         ifl:true=asymmetric
!
  integer p
  integer q
  integer r
  integer t

  real o(n,n),d(n,n),x(n,6),pijp(n,n),a(6,6)
  logical ifl
  logical*1 xceed(6,6)
!
  d(1:n,1:n) = 0.0

  do p=1,r
    do q=1,t
      xceed(p,q)=.false.
    end do
  end do

  write(iout,6)(q,q=1,t)
6     format(1x,'derivation of the generalized weight matrix.', &
     //' initial direction coefficients'/' dimension        direction'/ &
     3x,(7i10))

  do p=1,r
    write(iout,7)p,(a(p,q),q=1,t)
  end do

7     format(i6,(7f10.4))
  knt=1

  if ( ifull /= 1) then
    write(iout,44)
  else
    write(iout,45)
  end if

44    format(/' history of iterations'/'   iter   change   sstress   improve')
45    format(/' history of iterations'/'   iter   change   sstress', &
      '   improve',t42,'delta  q  p knt')
!
! calculate initial distances and sstress
!
  rsq=0.0
  call pddist(x,a,n,r,t,d)
  call pdstrs(o,d,n,ifl,oldfit,rsq)
  fitb4=oldfit
  stress=sqrt(oldfit/sumsq)
  write(iout,55)stress
55    format(5x,'0',3x,'initial',f10.5)
!
! perform als iterations
!
  kntcal=0
  do 50 iter=1,itmax
    big=0.0
    do 40 q=1,t
      do 30 p=1,r
        if ( xceed(p,q))go to 30
        call pdwght(a,x,d,o,n,r,pijp,p,q,ifl,delta)
        big=max (abs(delta),big)
        if ( abs(delta) < crit)xceed(p,q)=.true.
        if ( ifull /= 1)go to 30
        kntcal=kntcal+1
        call pdstrs(o,d,n,ifl,fit,rsq)
        stress=sqrt(fit/sumsq)
        diff=fitb4-fit
        fitb4=fit
    write(iout,9)iter,big,stress,diff,delta,q,p,kntcal
9       format(i6,4f10.5,3i3)
30      continue
40      continue
    call pdstrs(o,d,n,ifl,fit,rsq)
    diff=oldfit-fit
    oldfit=fit
    stress=sqrt(fit/sumsq)
    write(iout,4)iter,big,stress,diff
5       format(i6,(7f10.5))
4       format(i6,3f10.5)
    if ( big > crit)knt=0
  if ( big > crit)go to 50
    knt=knt+1
  if ( knt == 2)go to 60
    do 122 p=1,r
    do 122 q=1,t
122     xceed(p,q)=.false.
50      continue
60      rsq=-1.0
    call pdstrs(o,d,n,ifl,stresk,rsq)
    write(iout,11)stresk,rsq
11      format(/' fit measures (of distances to disparities)'/ &
    ' kruskal''s stress (1) =',f6.3/' squared correlation  =',f6.3)

  return
end
subroutine pdstrs ( o, d, n, ifl, stress, rsq )
!
!*******************************************************************************
!
!! PDSTRS calculates unnormalized stress.
!
!
!        (formula 1,note 1)
! written by Forrest W. Young, 11 dec 79
! final change 23may83 fwy
!
! o:dissimilarity data(or disparities)
! d:distances in principal direction space
! n:number of observations (stimuli)
! ifl:.true=asymmetric

   real o(n,n),d(n,n)
   logical ifl
!
! lower half matrix
!
   stress=0.0
   na=2
   if ( ifl)na=1
   if ( rsq /= -1.0)go to 20
   do 10 i=na,n
   nb=i-1
   if ( ifl)nb=n
   do 10 j=1,nb
   if ( i == j)go to 10
   if ( o(i,j) > 0.0)o(i,j)= sqrt( o(i,j))
   if ( o(i,j) < 0.0)o(i,j)=-sqrt(-o(i,j))
   if ( d(i,j) > 0.0)d(i,j)= sqrt( d(i,j))
   if ( d(i,j) < 0.0)d(i,j)=-sqrt(-d(i,j))
10     continue
20     continue

  do i=na,n
    nb=i-1
    if ( ifl)nb=n
    do j=1,nb
      if ( i /= j ) then
        stress=stress+(o(i,j)-d(i,j))**2
      end if
    end do
  end do

   if ( rsq /= -1.0)return
   denom=0.0
   do 35 i=na,n
   nb=i-1
   if ( ifl)nb=n
   do 35 j=1,nb
   if ( i == j)go to 35
   denom=denom+o(i,j)**2
35     continue
   stress=sqrt(stress/denom)
   nele=n*(n-1)/2
   if ( ifl)nele=n*(n-1)
   prodod=0.0
   sumo=0.0
   sumsqo=0.0
   sumd=0.0
   sumsqd=0.0
   do 40 i=na,n
   nb=i-1
   if ( ifl)nb=n
   do 40 j=1,nb
   if ( i == j)go to 40
   sumd=sumd+d(i,j)
   sumo=sumo+o(i,j)
   sumsqo=sumsqo+o(i,j)*o(i,j)
   sumsqd=sumsqd+d(i,j)*d(i,j)
   prodod=prodod+d(i,j)*o(i,j)
40     continue
   rsq=(prodod-sumd*sumo/nele)**2/ &
     ((sumsqd-sumd*sumd/nele)*(sumsqo-sumo*sumo/nele))

  return
end
subroutine pdtuck ( wk, x, uk, n, r )
!
!*******************************************************************************
!
!! PDTUCK obtains Tucker's oblique decomposition of WK.
!
!
  integer r
  real uk(n,6),x(n,6),wk(6,6)

  do i=1,r
    wk(i,i)=1.0/sqrt(wk(i,i))
  end do

  do i=2,r
    do j=1,i-1
      wk(i,j)=wk(i,i)*wk(i,j)*wk(j,j)
      wk(j,i)=wk(i,j)
    end do
  end do

  do j=1,r
    wk(j,j)=(1.0/wk(j,j))**2
    do i=1,n
      uk(i,j)=x(i,j)*wk(j,j)
    end do
  end do

  return
end
subroutine pdwght ( a, x, d, o, n, r, pijp, p, q, ifl, delta )
!
!*******************************************************************************
!
!! PDWGHT calculates a(p,q), the weight of variable p on principal direction q.
!
!
! written by Forrest W. Young and cynthia h. null, 10 dec 79
!
! a: principal direction weights
! x: multivariate data (common space)
! d: distances in the principal directions space
! o: similarity data

! n:number of observations(stimuli)
! r:number of variables(attributes)
! p:variable whose weight is being obtaineded
! q:principal directions whose weight is being obtained
! ifl:true=asymmetric
! delta:change in value of a(p,q) from entry to exit

  logical ifl
  integer r,p,q
  real o(n,n),d(n,n),pijp(n,n),x(n,6),a(6,6)
!
  c0=0.0
  c1=0.0
  c2=0.0
  c3=0.0
  do i=2,n
    do j=1,i-1
      pijp(j,i)=0.0
      pijp(i,j)= x(i,p)-x(j,p)
!
! compute sum in equation 16 of note 1
! (upper triangle of pijp is sum)
!
      do ia = 1, r
        if ( ia /= p ) then
          pijp(j,i)=pijp(j,i)+(x(i,ia)-x(j,ia))*a(ia,q)
        end if
      end do
!
! compute equation 16 note 1
!
      dpqk=2.0*pijp(i,j)*a(p,q)*pijp(j,i)+pijp(i,j)**2*a(p,q)**2
!
! compute equation 18 of note 1 as modified to permit
! asymmetric data (ifl=.true) (fixed mar 1982)
!
      rpqk=o(i,j)-d(i,j)+dpqk
      if ( ifl)rpqk=(rpqk+o(j,i)-d(i,j)+dpqk)*.5
!
! compute cubic equation coeficients as in eq.23 of note 1
!
      c0=c0+pijp(i,j)*rpqk*pijp(j,i)
      c1=c1+pijp(i,j)**2*(rpqk-2.0*pijp(j,i)**2)
      c2=c2+pijp(i,j)**3*pijp(j,i)
      c3=c3+pijp(i,j)**4

    end do
  end do

  c0=-c0/c3
  c1=-c1/c3
  c2=3.0*c2/c3
!
! save old a(p,q) in apq and solve for new value
!
  apq=a(p,q)
  call scube(c2,c1,c0,a(p,q))
  delta=a(p,q)-apq
!
! update distance using change in a(p,q)
! substituted into equation 16 of note 1.
!
  do i=2,n
    do j=1,i-1
      d(i,j)=d(i,j)+2.0*pijp(i,j)*pijp(j,i)*delta &
        -pijp(i,j)**2*(apq**2-a(p,q)**2)
    end do
  end do

  return
end
subroutine pdyung ( ck, wk, xk, x, n, r, t, roots )
!
!*******************************************************************************
!
!! PDYUNG obtains Young's principal directions solution.
!
!
!  Parameters:
!
!       at entry
!           ck  is the (nonorthogonal) optimal subspace coefficients.
!
!           x   is the group space.
!
!           n   is the number of stimuli.
!
!           r   is the number of group space dimensions.
!
!           t   is the number of principal directions
!
!       at return
!           ck  is the orthogonalized optimal subspace coefficients
!
!           wk  is the subject's weight matrix (the transformation
!               applied to the group space to obtain the
!               principal directions personal subspace).
!
  integer r,t
  real xk(n,6),x(n,6),ck(6,6),wk(6,6),roots(6)
!
!  compute weights (transformation) according to equation 4, note 1.
!
  do i=1,r
    do j=i,r
      wk(i,j)=0.0
      do k=1,t
        wk(i,j)=wk(i,j)+ck(i,k)*ck(j,k)
      end do
      if ( i /= j)wk(j,i)=wk(i,j)
    end do
  end do
!
!  Find principal direction coefficients according to equations 5 and 6, note 1.
!
  call eigk ( wk, roots, r, 6 )

  do i=1,r
    do j=1,t
      ck(i,j)=wk(i,j)*sqrt(roots(j))
    end do
  end do
!
!  Restore destroyed weights
!
  do i=1,r
    do j=i,r
      wk(i,j)=0.0
      do k=1,t
        wk(i,j)=wk(i,j)+ck(i,k)*ck(j,k)
      end do
      if ( i /= j)wk(j,i)=wk(i,j)
    end do
  end do
!
!  Compute principal direction subspace according to equation 8, note 1.
!
  do i=1,n
    do j=1,t
      xk(i,j)=0.0
      do k=1,r
        xk(i,j)=xk(i,j)+x(i,k)*ck(k,j)
      end do
    end do
  end do

  return
end
subroutine plotr ( x, y, xa, ya, xi, yi, npoi, out, id, long )
!
!*******************************************************************************
!
!! PLOTR is a plot routine for ALSCAL.
!
!
!  the original plot routine by f.w. young has been completely rewritten
!  by b.erichson / a.bischoff, jan. 89, for running alscal on pc.
!  plotr now serves as a calling routine for the new plot routine pplot.
!
  integer long
!
  integer out
  real x(long)
  real y(long)
!
!  page size: nlines x ncolum
!
  nlines = 60
  ncolum = 80

  iout = out
  ih = nlines
  iw = ncolum

  call pplot ( x, y, npoi, xi, xa, yi, ya, id, iout, ih, iw )

  write (iout,32)
   32 format (t115,'x')

  return
end
subroutine polyf ( disp, dist, n, obs, r, alph, ndeg, ncst, nab )
!
!*******************************************************************************
!
!! POLYF does polynomial fitting.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  dimension disp(1),dist(1),obs(1)
  dimension r(5,5),xy(5),alph(5)
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,mdeg,nt,nbs,nbnbns
!
  nor=ndeg
  lin=0
  if ( nor == 1.and.ncst == 1) lin=1
  if ( lin == 1) nor=2
  nor1=nor+1

  do i=1,nor1
    xy(i)=0.0
    do j=1,i
      r(i,j)=0.0
    end do
  end do

  r(1,1)=nab
  do 12 k=1,n
  if ( obs(k) == big) go to 12
  g=1.0
  xy(1)=xy(1)+dist( k)
  l=0
  do 13 i=2,nor1
  i1=i-1
  do 13 j=i1,i
  g=g*obs(k)
  l=l+1
  r(i,j)=r(i,j)+g
  if ( l > nor) go to 13
  xy(l+1)=xy(l+1)+g*dist( k)
   13 continue
   12 continue
  if ( nor1 < 3) go to 117
  do 15 i=3,nor1
  k=i-2
  do 15 j=1,k
   15 r(i,j)=r(i-1,j+1)
  117 continue
  do 16 i=1,nor1
  do 16 j=1,i
   16 r(j,i)=r(i,j)
  if ( ncst == 1) go to 20
  do 21 i=2,nor1
  i1=i-1
  xy(i1)=xy(i)
  do 21 j=2,nor1
  j1=j-1
   21 r(i1,j1)=r(i,j)
  nor1=nor
   20 continue
  call minv(r,nor1,5)

  do i=1,nor1
    alph(i)=0.0
    do j=1,nor1
      alph(i)=alph(i)+r(i,j)*xy(j)
    end do
  end do

  if ( lin == 1) return
  do 25 k=1,n
  if ( obs(k) == big) go to 23
  disp( k)=alph(1)
  if ( ncst == 0) disp( k)=alph(1)*obs(k)
  if ( nor1 < 2) go to 18
  g=1.0
  if ( ncst == 0) g=obs(k)
  do 19 i=2,nor1
  g=g*obs(k)
   19 disp( k)=disp( k)+g*alph(i)
  go to 18
   23 disp(k)=big
   18 continue

! the following line has been added july 1982.  before this change
! all numerical levels of measurement were unconstrained except
! linear interval, which was constrained to yield positive estimates.
! after this change all numerical transformations except linear interval
! use piece-wise linear regression.  nefatively extimated disparities
! are replaced with zero estimates. change made by Young and sarle.
! implemented by brooks.

  if ( disp(k) < 0.0)disp(k)=0.0
   25 continue
  return
end
subroutine pplot ( x, y, npoint, x1, x2, y1, y2, modus, iout, ih, iw )
!
!*******************************************************************************
!
!! PPLOT is a plot routine for ALSCAL.
!
!
!  Discussion:
!
!    The routine generates a printer plot of array -x- vs. array -y-.
!
!  Parameters:
!
!    x(npoint) : x-coordinates
!    y(npoint) : y-coordinates
!       x1, x2 : bounds for x-axis
!       y1, y2 : bounds for y-axis
!       x-bounds are generated if x1 = x2
!       y-bounds are generated if y1 = y2
!    modus > 0 : axes will be included
!          < 0 : no axes will be included
!          ñ 1 : points will be counted
!          ñ 2 : points will be identified
!         iout : output unit
!           ih : length of page
!           iw : width  of page
!  if vector x equals vector y then points will be plotted along
!  the horizontal axis (no axes will be plotted)
!  there can be one or two sets of points:
!   ncol = nrow = npoint  ->  one set :   x  y      (plot:  a,b,c,...)
!   ncol + nrow = npoint  ->  two sets:
!                     - column stimuli    x1 y1     (plot:  a,b,c,...)
!                                         -----
!                     - row stimuli       x2 y2     (plot:  a,b,c,...)
!
!  Define maximum size of plot field.
!
  integer, parameter :: maxx = 101
  integer, parameter :: maxy = 101
  integer, parameter :: mxy = maxx * maxy
  integer, parameter :: maxx2 = maxx + 14
!
!  define no. of intervals on x-axis and y-axis ( 6 <= nint <= 10 )
!
  integer, parameter :: nintx = 10
  integer, parameter :: ninty = 10

  logical axes
  logical onedim
  logical twoset

  character alfa(52)
  character ( len = 52 ) alfa1
  character alfal(26)
  character alfas(26)
  character, parameter :: bar = '-'
  character char1
  character, parameter :: cross = '0'
  character feld(maxx,maxy)
  character ( len = 13 ) fform
  character, parameter :: ibar = '|'
  character ( len = 10101 ) feld1
  character ( len = 101 ) feld2(maxy)
  character line(maxx2)
  character ( len = 115 ) line1
  character, parameter :: lio = '+'
  character, parameter :: liu = '+'
  character, parameter :: reo = '+'
  character, parameter :: reu = '+'
  character symb(10)
  character ( len = 10 ) symb1
  character, parameter :: til = '|'
  character, parameter :: tir = '|'
  character, parameter :: tiu = ':'
  real x(npoint)
  real xvalue(nintx+1)
  real y(npoint)
!
  equivalence (feld,feld1,feld2),(line,line1),(symb,symb1)
  equivalence(alfa,alfa1),(alfa1(1:26),alfal),(alfa1(27:52),alfas)

  common /cpplot/nrow

  data symb1/'x23456789m'/
  data alfa1/'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'/
!
  axes =.false.
  twoset =.false.
  if (nrow < npoint) then
    twoset =.true.
    ncol=npoint-nrow
  end if
!
!  plot design
!  -----------
!  iw x ih : size of page           (given from calling routine)
!  nx x ny : size of plot field
!    nintx : intervals on x-axis
!    ninty : intervals on y-axis
!      iix : spaces per interval on x-axis
!      iiy : spaces per interval on y-axis
!+----------------------------------------------------------------------------+
!| header                                                                     |
!| ------                                                                     |
!|<----14----><------------------------- nx x maxx -----------------------><2>|
!|           +----------------------------------------------------------------+  |
!|       1.0 |                              |    i                         |  |
!|       0.8 |                              |    i                         |  |
!|       0.6 |     plot field               |    i ny x maxy               |  |
!|       0.4 |                              |    i                         |  |
!|       0.2 |                              |    i                         |  |
!|       0.0 ---------------------------------------------------------------  |
!|      -0.2 |                              |    i                         |  |
!|      -0.4 |                              |    i                         |  |
!|      -0.6 |                              |    i                         |  |
!|      -0.8 |                              |    i                         |  |
!|      -1.0 |                              |    i                         |  |
!|           +----------------------------------------------------------------+  |
!|          -1.0  -0.8  -0.6  -0.4  -0.2   0.0   0.2   0.4   0.6   0.8   1.0  |
!+----------------------------------------------------------------------------+
!  size of plot field
!  ------------------
  nx = iw - 15
  ny = ih -  6
  nx = min ( maxx, nx )
  ny = min ( maxy, ny )
  iix = (nx-1)/nintx
  iiy = (ny-1)/ninty
  nx = iix * nintx + 1
  ny = iiy * ninty + 1
  midx = nx/2 + 1
  midy = ny/2 + 1
!
!  clear plot field
!
  feld1=' '
  line1=' '
!
!  check to see if one dimensional plot
!
  onedim=.false.
  do i=1,npoint
    if ( x(i) /= y(i))go to 110
  end do
  onedim=.true.
!
!  if desired, include axes on center of plot
!
  110 continue

  if (.not.onedim) then
    if (modus > 0) then
      axes=.true.
      do i=1,ny
        feld(midx,i)=ibar
      end do
      do i=1,nx
        feld(i,midy)=bar
      end do
      feld(midx,midy) = cross
    end if
  end if
!
!  determine minima and maxima, if necessary
!
!  x-axis
!  make range symmetric around zero
!
  xmin=x1
  ymin=y1
  xmax=x2
  ymax=y2

  if (xmax == xmin) then
    xmin = x(1)
    xmax = xmin
    do i = 2,npoint
      if (x(i) < xmin) xmin=x(i)
      if (x(i) > xmax) xmax=x(i)
    end do
    if (axes) then

      zmax = max (abs(xmin),abs(xmax))
      xmin = -zmax
      xmax =  zmax
    end if
  end if
!
!  y-axis
!  make range symmetric around zero
!
  if (ymax == ymin) then

    ymin = y(1)
    ymax = ymin
    do i = 2,npoint
      if (y(i) < ymin) ymin=y(i)
      if (y(i) > ymax) ymax=y(i)
    end do

    if ( axes ) then
      zmax = max ( abs ( ymin ), abs ( ymax ) )
      ymin = - zmax
      ymax =  zmax
    end if

  end if

  if (iabs(modus) == 1) then
  if (xmin <= 0.0) xmin=.00001
  if (ymin <= 0.0) ymin=.00001
  end if
!
!  Determine range and increment of both axes
!
  spanx = xmax-xmin
  spany = ymax-ymin
  delx  = spanx/(nx-1)
  dely  = spany/(ny-1)
!
!  Enter points into plot.
!
  do ii = 1, npoint
!
!  Locate point: (x,y) -> (j,i)
!
    if (onedim) then
      i=midy
    else
      i=(ymax-y(ii))/dely+1.5
      if (i > ny.or.i < 1) go to 200
    end if

    j=(x(ii)-xmin)/delx+1.5
    if (j > nx.or.j < 1) go to 200
!
!  identify point
!
    if ( iabs(modus) == 2) then

      if (twoset) then
        if ( ii <= ncol) then
          kk=mod(ii-1,26)+1
          feld(j,i)=alfas(kk)
        else
          kk=mod(ii-ncol-1,26)+1
          feld(j,i)=alfal(kk)
        end if
      else
        kk=mod(ii-1,52)+1
        feld(j,i)=alfa(kk)
      end if

    else
!
!  count point
!
      char1=feld(j,i)

      if ( char1 == ' ' .or. char1 == bar .or. char1 == ibar .or. &
           char1 == cross ) then
        feld(j,i)=symb(1)
      else
        do jj=1,10
          if (char1 == symb(jj)) then
            jj1=min ( jj+1, 10 )
            go to 202
          end if
        end do
  202       feld(j,i)=symb(jj1)
      end if
    end if

  200   continue

  end do
!
! print plot
! ----------
! top of frame
!
  write ( iout, * ) ' '
  write ( iout, * ) ' '
  write ( iout, * ) ' '

  do i = 14, nx+13
    line(i) = bar
  end do

  line(4)='Y'
  line(13)=lio
  line(nx+14)=reo
  line(midx+13)=tiu
  write(iout,300) line1(1:nx+14)
  300 format(1x,a)
!
!  cyclus over lines
!
  value = ymax+dely

  do i = 1, ny

    value=value-dely
    l=i+iiy-1
    if (l/iiy*iiy == l) then
      if (axes.and.i == midy) then
        write(iout,311) value,cross,feld2(i)(1:nx),cross
      else
        write(iout,311) value,til,feld2(i)(1:nx),tir
      end if
  311     format(' ',f11.1,' ',4a)
    else
      write(iout,'(13x,3a)') ibar,feld2(i)(1:nx),ibar
    end if

  end do
!
!  bottom of frame
!
  do j=1,nx
    l=j+iix-1
    if (l/iix*iix == l) line(j+13)=tiu
  end do

  line(4)=' '
  line(13)=liu
  line(nx+14)=reu
  line(midx+13)=cross
  write(iout,300) line1(1:nx+14)
!
!  labels for x-axis
!
  xvalue(1) = xmin
  do i=1,nintx
    xvalue(i+1)=xvalue(i)+iix*delx
  end do
!
!  format for labels (e.g.: '(10x,11f06.1)')
!
  iix1=14-iix+2
  write (fform,331) iix1,iix
  331 format ('(',i2,'x,11f',i2,'.1)')
  write(iout,fform) (xvalue(i),i=1,nintx+1)

  return
end
subroutine r_swap ( x, y )
!
!*******************************************************************************
!
!! R_SWAP switches two real values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, real X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  real x
  real y
  real z
!
  z = x
  x = y
  y = z

  return
end
subroutine prs ( ir, ibk, d, n, w, iz, l )
!
!*******************************************************************************
!
!! PRS prepares ties for continuous ordinal transformation.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     arrange tied observations in order of estimated values
!
  integer ir,ibk,iz
  dimension ir(1),ibk(1),d(1),w(1),iz(1)
!
  nb = 1
  if ( l /= 1 ) then
    nb=iz(l-1)
  end if

  ne = iz(l) - 2

  if ( ne < nb ) then
    return
  end if

  do jj = nb, ne, 2

    i=ibk(jj)
    if ( i > n) return
    j=ibk(jj+1)

    do i2=1,j
      k=i+i2-1
      w(i2)=d(ir(k))
    end do

    call shel9(w(1),ir(i),j)

  end do

  return
end
subroutine scube ( p, q, r, cube ) 
!
!*******************************************************************************
!
!! SCUBE solves a cubic equation.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  real, parameter :: r3 = 1.0 / 3.0
  real, parameter :: r27 = 3.703704E-02
!
  f(x)=(((x*0.25+p*r3)*x+q*0.5)*x+r)*x

  a = (3.0*q-p**2)*r3
  b = (p**2+p**2-9.0*q)*p*r27+r
  d = b**2*0.25+a**3*r27

  if ( d < 0.0 ) go to 1
  e = sqrt ( d ) - 0.5 * b
  cube = sign ( abs(e)**r3, e ) - p * r3
  e = - ( e + b )
  cube = sign(abs(e)**r3,e) + cube
  return

1 continue

  t3=acos(-b/(2.0*sqrt(-a**3*r27)))*r3
  e=2.0*sqrt(-a*r3)

  c1=cos(t3)
  c2=cos(t3+4.1887902)
  c3=cos(t3+2.0943951)

  cube= min (c1,c2,c3)*e-p*r3
  cmax= max (c1,c2,c3)*e-p*r3

  if ( f(cmax) < f(cube) ) then
    cube = cmax
  end if

  return
end
subroutine ses ( ir, ibk, d, n, iz, l )
!
!*******************************************************************************
!
!! SES prepares ties for discrete ordinal transformation.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     take the mean of tied observations
!
  integer ir,ibk,iz
  dimension ir(1),ibk(1),d(1),iz(1)

  nb=1
  if ( l /= 1 ) then
    nb=iz(l-1)
  end if

  ne=iz(l)-2

  if ( ne < nb ) then
    return
  end if

  do ii = nb, ne, 2

    i=ibk(ii)
    if ( i > n) return
    j=ibk(ii+1)
    jj=i+j-1

    s=0.0
    do k=i,jj
      s=s+d(ir(k))
    end do

    s=s/j
    do k=i,jj
      d(ir(k))=s
    end do

  end do

  return
end
subroutine shel9 ( a, c, nitem )
!
!*******************************************************************************
!
!! SHEL9 sorts data using Shell's sort.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     subroutine shel9   ricardo dobson
!     a is the key vector, c is the vector to sort on
!     nitem is the number of items in the two vectors
!     first vector argument must be real, second vector arg
!     must be integer
!     sort will be in ascending order
!
  real a(1)
  integer c(1)
  integer kk
  integer m
  integer nitem
!
  m = nitem
20    continue
  m = m / 2
  if ( m ) 30,40,30
30    k=nitem - m
  j=1
41    i=j
49    l=i+m
  if ( a(l)-a(i))50,60,60
50    b = a(i)
  a(i)= a(l)
  a(l)= b
  kk=c(i)
  c(i)=c(l)
  c(l)=kk
  i=i-m
  if (i-1)60,49,49
60    j=j+1
  if ( j-k)41,41,20
40    continue
  return
end
subroutine step0 ( nwords, eoj )
!
!*******************************************************************************
!
!! STEP0 reads and checks the problem parameters.
!
!
!  Discussion:
!
!    This routine reads and checks the problem parameters.
!    It then calculates the storage extent and structure for the problem.
!    The common /starts/ holds the starting positions of various
!    arrays used by the problem within a large block of storage
!    allocated by the main routine.  The naming convention used is
!    to prefix the array name by the letter j. ie. jxx gives the
!    starting position of the array xx  etc.
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  character ( len = 40 ) ccc
  character ( len = 80 ) fmt
  character ( len = 80 ) title
!
  common /cpplot/nrow
  common /pdcom/iflpds,ndir
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /block3/ title,fmt
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
  common /starts/jx,jwa,jwd,jxx,jds,jcfr,jcfl,jw,jws,jzz,jtr,jfk, &
         jcv,jcw,ju11,ju12,ju22,jub1,jub2, &
         jxn,jphsub,jphsti,jix,jiy,jiz,jndsr,jnad, &
         jftln,jdist,jpijp,jxk
  common /inicon/initx,initw,initws,initxc
  integer eoj

  double precision cut,stmin,stmind
  integer debug,icnstr,noulb
  common /prmblk/cut,stmin,debug,icnstr,noulb
!
!      fixed constants
!      ---------------
!
!  maxdim=6 defines the maximum number of dimensions.
!      this value may not be changed.
!  big=9.0e20 defines the largest nominal or ordinal datum allowed.
!      this may be changed.
!  epsid=.001 defines the default iteration convergence criterion.
!      this may be changed.
!  maxitd=30 defines the default maximum number of iterations.
!      this may be changed.
!  stmind=.005 defines the default minimun s-stress.
!      this may be changed.
!
  parameter ( maxdim = 6 )
  parameter ( epsid = 0.001 )
  parameter ( maxitd = 30 )
  parameter ( stmind = 0.005 )

!-pc  data d,f,dyn,fix/'d','f','dynamic ','  fixed '/
  big=9.0e20

  eoj = 0

    2 nadct=0
!
!  CARD 1: job start/end
!
    1 continue

  read ( in, '(a)', end = 905 ) title

  if ( title == 'end' .or. title == ' end' .or. &
       title == 'END' .or. title == ' END' ) then

    write ( lout, * ) ' '
    write ( lout, * ) 'STEP0 - Note:'
    write ( lout, * ) '  END card encountered in input file.'

    write ( *, * ) ' '
    write ( *, * ) 'STEP0 - Note:'
    write ( *, * ) '  END card encountered in input file.'

    eoj = 1
    return

  end if
!
!     data specifications
!
!     nrow     number of row stimuli (no maximum)
!     ncol     number of column stimuli (no maximum)
!     ns       number of matrices (no maximum)
!     ndtyp    measurement level
!                   =1   ratio (polynomial w/o intercept)
!                   =2   interval (polynomial with intercept)
!                   =3   ordinal(default)
!                   =4   nominal
!     nsim     data form
!                   =0   symmetric-dissimilarity(default)
!                   =1   symmetric-similarity
!                   =2   asymmetric-dissimilarty
!                   =3   asymmetric-similarity
!                   =4   rectangular-dissimilarity
!                   =5   rectangular-similarity
!     nps      measurement process (only when ndtyp=3)
!                   =1   discrete(default)
!                   =2   continuous
!     nwc      measurement conditionality
!                   =1   matrix conditional(default)
!                   =2   row conditional
!                   =3   unconditional
!     ndeg     degrees in polynomial (when ndtyp=1 or 2)
!     ndmx     maximum number of ordinal observation categories.
!                   (default= total number of cells, or 1000,
!                             whichever is smaller)
!     cut      cutoff for missing data (default 0.0)
!
  read ( in, 101, err = 905, end = 905 ) nrow, ncol, ns, ndtyp, nsim, nps, &
    nwc, ndeg, ndmx, cut

  101 format(9i4,f8.4)
  if ( nrow < 3)go to 903
  if ( ns <= 0)go to 904

  call page ( lout )

  write ( lout, * ) 'ALSCAL'
  write ( lout, * ) '  Alternating Least Squares Scaling'
  write ( lout, * ) ' '
  write ( lout, * ) '  Yoshio Takane,'
  write ( lout, * ) '  Forrest W. Young and'
  write ( lout, * ) '  Rostyslaw Lewyckyj.'
  write ( lout, * ) ' '
  write ( lout, * ) '  Adaptation by B. Erichson and A. Bischoff'
  write ( lout, * ) ' '
  write ( lout, * ) '  Psychometric Laboratory'
  write ( lout, * ) '  The University of North Carolina'
  write ( lout, * ) '  Chapel Hill, NC  27514'
  write ( lout, * ) ' '
  write ( lout, * ) '  Copyright 1977'
  write ( lout, * ) '  F. W. Young, Y. Takane and R. J. Lewyckyj'
  write ( lout, * ) ' '
  write ( lout, * ) ' '

  nb=nrow
  if ( nsim > 3)nb=nb+ncol
  nbs=nb**2
  nbnbns=nb**2*ns
  nc=nb*(nb-1)
  nc2=nc/2
  ndx=min ( maxdim, nb-2 )

  if ( ns == 1) then
    ccc='matrix'
  else
    ccc='matrices'
  end if

  write(lout,200)title,nrow,ncol,ns,ccc
  200 format(//' job title: ',a80// ' data specifications-'/ &
    /'  nrow - number of row stimuli',t50,i5,t58,'row stimuli'/ &
    '  ncol - number of column stimuli',t50,i5,t58,'column stimuli'/ &
    '  ns   - number of matrices',t50,i5,3x,a8)
  ncst=1

  if ( ndtyp == 1) then
    ccc='ratio'
    ncst=0
  else if ( ndtyp == 2) then
    ccc='interval'
  else if ( ndtyp == 4) then
    ccc='nominal'
  else
    ccc='ordinal'
    ndtyp=3
  end if

  write(lout,208) ndtyp,ccc
  208 format('  ndtyp- measurement level',t50,i5,' = ',a8)

  if ( nsim < 0.or.nsim > 5) then
    nsim=0
  end if

  if ( nsim > 1) then
    nc2=nbs
  end if

  nt=nc2*ns

  if ( nsim == 0) then
    ccc='Symmetric-dissimilarity'
  else if ( nsim == 1) then
    ccc='Symmetric-similarity'
  else if ( nsim == 2) then
    ccc='Asymmetric-dissimilarity'
  else if ( nsim == 3) then
    ccc='Asymmetric-similarity'
  else if ( nsim == 4) then
    ccc='Rectangular-dissimilarity'
  else if ( nsim == 5) then
    ccc='Rectangular-similarity'
  end if

  write(lout,210)nsim,ccc
  210 format('  nsim - data type',t50,i5,' = ',a25)

  if ( nps == 2) then
    ccc='Continuous (untie)'
  else
    ccc='Discrete (tie)'
    nps=1
  end if

  write(lout,209) nps,ccc
  209 format('  nps  - measurement process',t50,i5,' = ',a20)

  if ( nps == 2.and.ndtyp == 4) then
    write(lout,9910)
    write(lout,228)
    write(*,9910)
    write(*,228)
    call hitr
    nps=1
  end if

  228 format(8x,'the program will continue with the discrete' &
    ,' measurement process (nps=1)')

  if ( nwc == 2) then
    ccc='Rowconditional'
  else if ( nwc == 3) then
    ccc='Unconditional'
    nwc=0
  else
    ccc='Matrix-conditional'
    nwc=1
  end if

  write(lout,261) nwc,ccc
  261 format('  nwc  - measurement conditionality',t50,i5,' = ',a20)

!-----the following two statements added 8/4/82
  write(lout,203)cut
  203 format('  cut  - data cutoff',t50,f13.7)

  if ( ndeg < 1) then
    ndeg=1
  end if

  if ( ndtyp <= 2) then
  if ( ndeg <= 1) go to 813
  if ( ndeg == 2) then
    ccc='Quadratic'
  else if ( ndeg == 3) then
    ccc='Cubic'
  else
    ccc='Quartic'
    ndeg=4
  end if

  write(lout,811)ndeg,ccc
  end if
  811 format('  ndeg - degrees in polynomial',t50,i5,' = ',a10)
  813 if ( nrow /= ncol.and.nsim < 4)go to 902
!
!     analysis specifications
!     -----------------------
!
!     nwe      model type
!                   =0   simple Euclidean model (default)
!                   =1   individual differences model
!                   =2   multiprocess asymmetric model
!                   =3   multiprocess asymmetric individual differences
!                        model
!                   =4   generalized Euclidean model
!     ndim     number of dimensions in the solution  (maximum)
!     ndmn     number of dimensions in the solution  (minimum)
!     nnc=1    if weights not constrained to be positive
!              (default is nonnegativity restrictions)
!     maxit    maximum number of iterations (default is 30)
!     epsi     convergence criterion (default is 0.001)
!     stmin    minimum stress (default is .005)
!     ndir     number of gemscal directions


!-----stmin added to following statement 7/7/82
!-----ndir added to following statement for gemscal 8/5/82
  read(in,102)nwe,ndim,ndmn,nnc,maxit,epsi,stmin,ndir
  102 format(5i4,2f8.0,i4)

  write(lout,216)
  216 format(/'Analysis specifications-'/)

  if ( nwe == 1) then
    ccc='Individual differences (INDSCAL) model'
  else if ( nwe == 2) then
    ccc='Asymmetric Euclidean (ASYMSCAL) model'
  else if ( nwe == 3) then
    ccc='Asymmetric INDSCAL (ASYNDSCAL) model'
  else if ( nwe == 4) then
    ccc='Generalized Euclidean (GEMSCAL) model'
  else
    ccc='Simple Euclidean model (Default)'
    nwe=0
  end if

  write(lout,218) nwe,ccc
  218 format('  nwe  - model type',t50,i5,'   ',a40)

  iflpds=0
  if ( nwe == 4) then
    iflpds=1
    nwe=1
  end if

  if ( nwe == 1.and.ns == 1) then
    write(lout,9910)
    write(lout,333)
    write(*,9910)
    write(*,333)
    call hitr
  end if

  333 format(8x,'the individual differences or generalized Euclidean', &
     ' model',8x,'can not be used with only one matrix.'/ &
     8x,'analysis continues with Euclidean model (nwe=0)')

  if ( nwe == 1.and.ns == 1) nwe=0

  if ( nwe == 3.and.ns == 1) then
    write(lout,9910)
    write(lout,9238)
    write(*,9910)
    write(*,9238)
    call hitr
  end if

 9238 format(8x,'the individual differences asymmetric model cannot', &
     ' be used with only one matrix.'/ &
     8x,'analysis continues with the asymmetric model (nwe=2)')
  if ( nwe == 3.and.ns == 1) nwe=2
  if ( nwc /= 0.or.(nwe /= 0.and.nwe /= 2).or.ns /= 1)go to 17
  nwc=1
  write(lout,9910)
  write(lout,219)
  write(*,9910)
  write(*,219)
  call hitr
  219 format(8x,'program will continue treating the data as matrix', &
     ' conditional (nwc=1)')
   17 continue
  if ( nwe <= 1.or.nsim > 1)go to 20
  write(lout,9910)
  write(lout,21)
  write(*,9910)
  write(*,21)
   21 format(8x,'the asymmetric model requires asymmetric data.', &
     /8x,' the model has been changed as follows')
  nwe=nwe-2

  if ( nwe == 0) then
    ccc='simple Euclidean model (default)'
  else
    ccc='individual differences (indscal) model'
  end if

  write(lout,218)nwe,ccc
  write(*,218)nwe,ccc
  call hitr

   20 if ( nwc /= 2.or.nsim >= 2)go to 22
  write(lout,9910)
  write(lout,23)
  write(*,9910)
  write(*,23)
  call hitr
   23 format(8x,'row conditional data are not permitted with symmetric data.', &
     /8x,'the data will be treated as matrix conditional. (nwc=1)')
  nwc=1

   22 continue

  if ( nwe > 1.and.nwc == 2) then
    write(lout,9910)
    write(lout,9901)
    write(*,9910)
    write(*,9901)
    call hitr
    nwc=1
  end if

 9901 format(8x,'row conditional data are not permitted with the', &
     ' asymmetric models.'/ &
     8x,'analysis continues with matrix  conditional data (nwc=1)')
  if ( nwe > 1.and.ndim > 5) ndim=5
  if ( ndim == 0) ndim=2

  if ( ndim > ndx) then
    write(lout,9910)
    write(lout,4003) ndx
    write(*,9910)
    write(*,4003) ndx
    call hitr
  end if

 4003 format(8x,'the maximum dimensionality may not exceed',i2)

  if ( ndim > ndx) ndim=ndx
  ndx=ndim
  ndxs=ndx**2
  ndxp=ndx+1
  if ( ndmn <= 0.or.ndmn > ndim)ndmn=ndim

  if ( ndmn == 1.and.nwe > 0) then
    write(lout,9910)
    write(lout,9902)
    write(*,9910)
    write(*,9902)
    call hitr
  end if

  if ( ndmn == 1.and.nwe > 0.and.ndim == 1)ndim=2
  if ( ndmn == 1.and.nwe > 0)ndmn=2
 9902 format(8x,'one-dimensional weighted models not permitted.'/8x, &
    'analysis continues without a one-dimensional solution.')
  write(lout,217) ndim,ndmn
  217 format('  ndim - number of dimensions (maximum)',t50,i5, &
     t58,'dimensions (maximum)'/ &
     '  ndmn - number of dimensions (minimum)',t50,i5, &
     t58,'dimensions (minimum)')
  nd=ndim-ndmn+1
  nsc=1
  if ( nwc == 0)nsc=0
!-----following two statements added for gemscal 8/5/82
  if ( iflpds == 1)write(lout,8366)ndir
 8366 format('  ndir - number of gemscal directions ',t50,i5)

  if ( nnc == 1) then
    ccc='negative weights permitted'
  else
    ccc='negative weights not permitted'
    nnc=0
  end if

  write(lout,4532)nnc,ccc
 4532 format('  nnc  - negative weights permitted ',t50,i5,' = ',a30)
!
!     i/o options
!     -----------
!
!     ndt      print input data
!                  =0 no
!                  =1 yes
!     npt      plot results
!                  =0 no
!                  =1 plot spaces and overall fit
!                  =2 also plot transformation and fit for every
!                     partition (can be very many pages of output)
!     nph      punch results
!                  =0 do not punch
!                  =1 punch derived configuration
!                  =2 punch initial and derived configuration
!     indata   data input unit
!                  =0   data read from cards
!                  =n   data read from unit n
!     initx    initial stimulus coordinates
!                  =0   compute
!                  =1   compute and print
!                  =2   read and print
!                  =3   read, print and fix
!     initxc   initial column stimulus coordinates
!                  =0   compute
!                  =1   compute and print
!                  =2   read and print
!                  =3   read, print and fix
!     initw    initial subject weights
!                  =0   compute
!                  =1   compute and print
!                  =2   read and print
!                  =3   read, print and fix
!     initws   initial stimulus weights
!                  =0   compute
!                  =1   compute and print
!                  =2   read and print
!                  =3   read, print and fix
!     noulb    upper and lower bounds to estimate missing data
!                  =0   yes (default)
!                  =1   no  (alscal 4 method)
!     icnstr   constrain missing data (not fully implemented)
!     debug    full debugging output

!-----noulb, icnstr, and debug added to following statement 7/7/82
  read(in,103)ndt,npt,nph,indata,initx,initxc,initw,initws,noulb,icnstr,debug
  103 format(12i4)
  icnstr=0
  write(lout,204)
  204 format(/' i/o options-'/)
  if ( debug /= 0)write(lout,4387)
 4387 format(' debugging is turned on')

  if ( ndt == 0) then
  ccc='do not print'
  else
  ccc='do print'
  end if
  write(lout,205) ndt,ccc
  205 format('  ndt  - print data, distances and disparities',t50,i5,' = ',a15)

  if ( npt == 0) then
    ccc='do not plot'
  else
    ccc='do plot'
  end if

  write(lout,206) npt,ccc
  206 format('  npt  - plot results ',t50,i5,' = ',a15)

  if ( nph == 0) then
    ccc='do not punch'
  else if ( nph == 1)then
    ccc='punch derived configuration'
  else
    ccc='punch initial and derived configurations'
  end if

  write(lout,207) nph,ccc
  207 format('  nph  - punch results ',t50,i5,' = ',a40)

  if ( nph  >=  1 ) then

    open ( unit = nplt, file = 'nplt' )

  end if

  if ( indata == 0)indata=in
  if ( indata == in) then
    ccc='read data from cards'
  else
    ccc='read data from disk or tape'
  end if
  write(lout,243)indata,ccc
  243 format(' indata- data input unit number',t50,i5,' = ',a30)

  if ( initx == 0) then
    ccc='compute'
  else if ( initx == 1) then
    ccc='compute and print'
  else if ( initx == 2) then
    ccc='read and print'
  else
    ccc='read, print and fix'
  end if

  write(lout,211)initx,ccc
  211 format(' initx - initial stimulus coordinates',t50,i5,' = ',a20)

  if ( initxc == 0) then
    ccc='compute'
  else if ( initxc == 1) then
    ccc='compute and print'
  else if ( initxc == 2) then
    ccc='read and print'
  else
    ccc='read, print and fix'
  end if

  write(lout,244)initxc,ccc
  244 format(' initxc- initial column stimulus coordinates',t50,i5,' = ',a20)

  if ( initw == 0) then
    ccc='compute'
  else if ( initw == 1) then
    ccc='compute and print'
  else if ( initw == 2) then
    ccc='read and print'
  else
    ccc='read, print and fix'
  end if
  write(lout,222)initw,ccc
  222 format(' initw - initial subject weights',t50,i5,' = ',a20)

  if ( initws == 0) then
    ccc='compute'
  else if ( initws == 1) then
    ccc='compute and print'
  else if ( initws == 2) then
    ccc='read and print'
  else
    ccc='read, print and fix'
  end if
  write(lout,232)initws,ccc
  232 format(' initws- initial stimulus weights',t50,i5,' = ',a20)

!-----algorithmic options section comprised of moved and new statements
  write(lout,4902)
 4902 format(////' algorithmic options-'/)
  if ( maxit == 0)maxit=maxitd
  write(lout,5002)maxit
 5002 format('  maxit- maximum number of iterations',t50,i5,t58, &
    'iterations (maximum)')
  if ( epsi <= 0)epsi=epsid
  write(lout,5001) epsi
 5001 format('  epsi - convergence criterion',t45,f10.7, &
    ' = minimum sstress improvement')
  if ( stmin <= 0.0)stmin=stmind
  write(lout,5003)stmin
 5003 format('  stmin- minimum sstress',t45,f10.7,' = minimum sstress cutoff')
  if ( noulb == 1)write(lout,5981)noulb
 5981 format('  noulb- initial missing data estimates',t50,i5,' = means')
  if ( noulb == 0)write(lout,5982)noulb
 5982 format('  noulb- initial missing data estimates',t50,i5,' = ulbounds')
  if ( ndmx <= 0)ndmx=min ( nt, 1000 )
  if ( ndtyp /= 3)ndmx=0
  if ( ndtyp == 3)write(lout,812)ndmx
  812 format('  ndmx - number of cells for tied observations',t50,i5,t58,'cells')
  if ( icnstr == 1)write(lout,5983)
 5983 format(' unfolding analysis is constrained')

!
!     calculation of storage requirements for arrays
!
!  in this section we calculate the number of words of storage
!  required for various arrays used by the problem. this is done
!  by arranging the arrays end to end and adding their lengths
!  together.  In this process at each step the starting place of
!  the next array is calculated by adding the length of the current
!  array to its starting point. thus the starting point of the last
!  array plus its length gives the size of the storage block needed.
!
!  As a convention the starting point of an array is given by
!  the name of the array prefixed by j . thus
!      jnext=jcurr + lcurr
!  where lcurr is usually an expression in the parameters of the pro-
!  blem. if more than one array is to share the same storage then
!      jnext=jcurr + (maximum of the lengths of arrays sharing sto-
!   rage with curr) so as not to overwrite the next array.
!     the starting points of the arrays are stored in the common
!  block /starts/.
!     the final calculated number of words of storage required is saved
!  in nwords.

!----------------------------------------------------------------------

  jx = 1
!
!   x is used as a 2*nb*ns array in inner.
!
  jwa = jx + max ( nt, 2*nb*ns )
  jwd=jwa+nbnbns
  jxx=jwd+max ( nt, ndxs*nb, nbs*ndx, 2778 )
!   xx is passed to am in init and used there as an ns*ns array
  jds=jxx+max ( nbs, ns*ns )
  jcfr=jds+nbs
  lcfr=nb*ndxp
  jcfl=jcfr+lcfr
!
!  cfl gets passed to xeq in init and from there in two pieces
!  to u and v in cjeig. for this call the length of u+v is 4*ns.
!  otherwise cfl is used as an nb*ndxp array.
!
  jw=jcfl+max ( 4*ns, lcfr )
  jws=jw+ns*ndxp
  jzz=jws+lcfr
!
!  zz gets passed from inswm to cjeig as scratch space of length 2*nb.
!  otherwise zz is used as a square (ndim*(ndim+1))/2 array.
  ju11=jzz+max ( 2*nb, ((ndim+1)*ndim/2)**2 )
  ju12=ju11+ndxs
  ju22=ju12+ndxs
  jub1=ju22+ndxs
  jub2=jub1+ndim
  jxn=jub2+ndim
  jtr=jxn+ndim
!  tr is used as scratch space in several subroutines
  jfk=jtr+max ( nb, ns )
  jcv=jfk+ndxp**2
  lcv=((ndim+1)*ndim)/2
  jcw=jcv+lcv
!
!       lcw=lcv
!
  jphsub=jcw+lcv
  jphsti=jphsub+ns
  jndsr=jphsti+nb

  if ( nwc == 1 .and. ndtyp == 4 ) then
    jnad = jndsr + ns
  else if ( nwc == 2 .and. ndtyp == 4 ) then
    jnad = jndsr + ns * nb
  else
    jnad = jndsr
  end if

  if ( nwc == 1 ) then
    jix = jnad + ns
  else if ( nwc == 2 ) then
    jix = jnad + ns * nb
  else
    jix = jnad
  end if

  jiy=jix+(nt+1)
  jiz=jiy+(ndmx+1)
  jftln=jiz+1
  if ( ndtyp /= 3)go to 5100
  if ( nwc == 1)jftln=jiz+(ns+1)
  if ( nwc == 2)jftln=jiz+(nb*ns+1)
!-----space for ftln not being used in version 4.04 and up
 5100 jdist=jftln+2
  if ( ndtyp /= 2.or.ndeg /= 1.or.nwc == 0)go to 510
  if ( nwc == 1)jdist=jftln+ns*2
  if ( nwc == 2)jdist=jftln+ns*nb*2

  510 if ( iflpds /= 1)go to 7381
  jpijp=jdist+nbs
  jxk=jpijp+nbs
  nwords=jxk+nb*6
  return

 7381 jpijp=jdist
  jxk=jdist
  nwords=jdist

  return

  902 write(lout,9900)
  write(*,9900)
  write(lout,9918)
  write(*,9918)
 9918 format (8x,'number of rows must equal number of columns for' // &
      ' non-rectangular data.')

  eoj=1
  return
  903 write(lout,9900)
  write(*,9900)
  write(lout,9926)
  write(*,9926)
 9926 format(8x,' number of stimuli less than 3'/)
  eoj=1
  return
  904 write(lout,9900)
  write(*,9900)
  write(lout,9927)
  write(*,9927)
 9927 format(8x,' total number of matrices less than 1'/)
  eoj=1
  return

  905 write ( lout, * ) ' '
  write ( lout, * ) 'ALSCAL message: All cards and problems read.'
  write ( lout, * ) '  Normal end of job.'
  eoj=1
  return
 9900 format(/' alscal fatal error:  computations terminated.')
 9910 format(/' alscal warning:  inconsistent control parameters.')
end
subroutine step1 ( ix, x, wa, xx, ds, nad, ier )
!
!*******************************************************************************
!
!! STEP1 is the data input routine.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  integer nb
  integer ns
  integer nt
!
  double precision ax
  double precision cut
  real ds(nb,nb)
  character ( len = 80 ) fmt
  integer ix(nt)
  integer nad(ns,nb)
  double precision stmin
  character ( len = 80 ) title
  real wa(nb,nb,ns)
  real x(nt)
  real xx(nb,nb)
!
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /block3/ title,fmt
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
!
  integer debug,icnstr,noulb
  common /prmblk/cut,stmin,debug,icnstr,noulb
!
!  x contains the active data ( with the missing data flagged)
!  in standard order with similarities converted to dissimilarities
!  and scaled to known units.
!  wa is the same as x except
!    1) has means in place of missing data
!    2) has additive constant estimated
!    3) is always symmetric
!  ds contains the sum,over subjects, of information in x
!  ix contains the missing data pattern 1=missing, 0=active
!  xx is used as temporary space to construct each matrix of
!        data to be used by the initialization routines.
!        xx is processed so that:
!           1) missing data estimated by mean upper-lower bound
!           2) asymmetric data are symmetrized
!           3) similarities are made dissimilarities
!           4) conditional data are normalized
!  wa contains the same information as xx, except
!           1) for all matrices
!           2) additive constant is estimated
!           3) result is squared
!
!  certain important data flags and their values are:
!      nsim: 0,1=sym 2,3=asym 4,5=rect even=dissim odd=similarities
!      nwc:  0=unconditional 1=matcon 2=rowcon
!
  read ( in, '(a80)' ) fmt
!   fmt is the variable format for the input data
  write(lout,211)fmt
  211 format(/' input data format-'//5x,a80)
  if ( ndt == 1) then
  call page ( lout )
  write(lout,212)
  212   format(//' input data')
  end if
  nab=0
  nx=2
  if ( nsim > 1) nx=1
  i1=nb
  do 4001 i=nx,nb
  if ( nsim <= 1)i1=i-1
  do 4001 j=1,i1
 4001 ds(i,j)=0.0
  nn=0
  ier=0
  rewind ndpp

  do 10 l = 1, ns

    if ( debug /= 0 ) then
      write(lout,4040)
    end if

 4040 format(/' raw data as read from data set')
    ntemp2=0
    xx(1,1)=0.0

    do i = nx, nb

      if ( nsim <= 1)i1=i-1
      if ( nsim > 3.and.i <= ncol)go to 7
      j1=i1
      if ( nsim > 3)j1=ncol
      read(indata,fmt)(xx(i,j),j=1,j1)
      if ( debug /= 0)write(lout,40)(xx(i,j),j=1,j1)
   40     format(6x,10g12.4)
      if ( nsim < 4)go to 9
      j1=j1+1
      do j = j1, i1
        xx(i,j) = cut - 1.0
      end do
      go to 9
    7     do j=1,nb
        xx(i,j)=cut-1.0
      end do
    9     write(ndpp)(xx(i,j),j=1,i1)

      if ( nsim <= 1 ) then

        do j = 1, i1
          xx(j,i) = xx(i,j)
        end do

      end if

    end do
!
!-----the following ten lines added 7/9/82
!*************the sas version reads in a transposed matrix and thus must
!     retranspose it.  this is not the case with the stand alone
!     version, so this section of code is not needed.*****

!     if the matrix is not symmetric, transpose it

!     if ( nsim < 2)go to 21
!     do 20 i=2,nb
!     do 20 j=1,i-1
!     temp=xx(i,j)
!     xx(i,j)=xx(j,i)
!  20 xx(j,i)=temp
!
   21   if ( ndt /= 1) go to 13
    write(lout,215) l
  215   format(/' matrix',i4/)
    if ( nsim <= 1)call outs(x,nb,1,-1,xx)
    if ( nsim > 1)call outa(xx,nb,-1)
13      continue
!
!  Check for missing data and replace by means.
!
    do i=1,nb
      xx(i,i)=0.0
    end do

    fmr=0.0
    n=nn
    iflxer=0

    if ( debug  >  0 ) then
      write ( lout, * ) 'Missing data pattern.'
    end if

    do 2337 i=nx,nb
    ntemp1=0
    t=0.0
    k=n
    if ( nsim <= 1)i1=i-1
    do 337 j=1,i1
    n=n+1
!---the following statement is replaced by the next to fix data cutoff
!       if ( xx(i,j) > 0.0.or.i == j)go to 338
    if ( xx(i,j) >= cut.or.i == j)go to 338
    ix(n)=1
!-----the following three statements added 7/9/82
    if ( icnstr /= 1.or.nsim < 4)go to 337
    if ( i > ncol.and.j <= ncol)go to 337
    ix(n)=0
    go to 337
  338   ntemp1=ntemp1+1
    t=t+xx(i,j)
    ix(n)=0
  337   continue
    if ( ntemp1 == 1.and.nwc == 2) iflxer=1
    if ( ntemp1 == 0.and.nwc == 2)write(lout,9001)l,i
 9001   format(/' alscal warning: all data missing for matrix',i4, &
        'stimulus',i4)
    ntemp2=ntemp2+ntemp1
    fmr=fmr+t
!-----the following two statements replaced by the next three 7/9/82
!       if ( nwc /= 2.or.iflxer == 1)go to 2337
!       nad(l,i)=ntemp1
    if ( nwc == 2)nad(l,i)=ntemp1
    if ( nwc /= 2)go to 2337
    if ( nwc == 2.and.ntemp1 == 1)go to 2337
    t=t/(ntemp1-1)
    do j=1,i1
      k=k+1
      if ( ix(k) /= 0)xx(i,j)=t
    end do

 2337   continue
!-----the following two statements added 7/9/82
    if ( debug > 0)write(lout,4789)(ix(iiiiii),iiiiii=1,n)
 4789   format(10(1x,10i1))
    nmiss=nc2-ntemp2
    if ( nsim > 3)nmiss=nmiss-ncol*nb-(nb-ncol)**2+nb

    if ( nmiss  >  0 ) then

      write ( lout, * ) ' '
      write ( lout, * ) 'ALSCAL message:'
      write ( lout, * ) '  For matrix ', l
      write ( lout, * ) '  Number of missing observations: ', nmiss

      write ( *, * ) ' '
      write ( *, * ) 'ALSCAL message:'
      write ( *, * ) '  For matrix ', l
      write ( *, * ) '  Number of missing observations: ', nmiss

      call hitr

    end if

    if ( ntemp2 == 0)go to 901
    if ( ier /= 0)go to 10
    nab=nab+ntemp2
!-----the following two statements added 7/9/82
    if ( nsim <= 1)fmr=fmr/ntemp2
    if ( nsim > 1)fmr=fmr/(ntemp2-nb)
    if ( nwc == 2.and.iflxer == 0)go to 341
!-----the following statement added 7/9/82
    if ( nwc == 2)go to 342
    if ( nwc == 1)nad(l,1)=ntemp2
    if ( nmiss == 0.and.nsim < 4)go to 341
    n=nn
    do 339 i=nx,nb
    if ( nsim <= 1)i1=i-1
    do 339 j=1,i1
    n=n+1
    if ( ix(n) /= 0)xx(i,j)=fmr
!-----the following statement added 7/9/82
    if ( ix(n) /= 0.and.nsim <= 1)xx(j,i)=fmr
  339   continue
!-----the following seven statements added 7/9/82
    go to 341
  342   n=nn

    do i=1,nb
      do j=1,nb
        n=n+1
        if ( nad(l,i) > 1)go to 343
        if ( ix(n) /= 0)xx(i,j)=fmr
  343       continue
      end do
    end do
!
!     if data are similarity (nsim=1 or 3) convert them
!     into dissimilarity
!
  341   if ( nsim/2*2 == nsim) go to 114
    if ( nwc == 2) go to 2005
    ax=xx(nb,1)

    do i=nx,nb
      if ( nsim <= 1)i1=i-1
      do j=1,i1
        if ( xx(i,j) > ax) ax=xx(i,j)
      end do
    end do

    do i=nx,nb
      if ( nsim <= 1)i1=i-1
      do j=1,i1
        xx(i,j)=ax-xx(i,j)
      end do
      xx(i,i)=0.0
    end do

    xx(1,1)=0.0
    go to 114
 2005   do 2006 i=1,nb
    ax=xx(i,1)
    do 2007 j=1,nb
 2007   if ( xx(i,j) > ax) ax=xx(i,j)
    do 2008 j=1,nb
 2008   xx(i,j)=ax-xx(i,j)
    xx(i,i)=0.0
 2006   continue
!
!     initial scaling within subject
!
!-----the following five statments were added 7/8/82
  114   if ( debug == 0) go to 3200
    write(lout,3399)
 3399   format(/' data with missing estimated by row or matrix mean')
    do 3388 i=1,nb
 3388   write(lout,6789)(xx(i,j),j=1,nb)
 3200   ax=1.0
    if ( nsc == 0) go to 5091

    ax=0.0
    do i=nx,nb
      if ( nsim <= 1)i1=i-1
      do j=1,i1
        ax=ax+xx(i,j)**2
      end do
    end do

    ax=dsqrt(ax/nc2)

    do i=nx,nb
      if ( nsim <= 1)i1=i-1
      do j=1,i1
        xx(i,j)=xx(i,j)/ax
      end do
    end do

    if ( nsim > 1)go to 5091

    do i=nx,nb
      i1=i-1
      do j=1,i1
        xx(j,i)=xx(i,j)
      end do
    end do
!
!     store original data
!
 5091   do 4002 i=nx,nb
    if ( nsim <= 1)i1=i-1
    do 4002 j=1,i1
 4002   ds(i,j)=ds(i,j)+xx(i,j)
!-----the following five statements were added 7/8/82
    if ( debug == 0) go to 3210
    write(lout,3377)
 3377   format(/' normalized data and means')
    do 3366 i=1,nb
 3366   write(lout,6789)(xx(i,j),j=1,nb)
 3210   n=nn

    do i=nx,nb
      if ( nsim <= 1)i1=i-1
      do j=1,i1
        n=n+1
        x(n)=xx(i,j)
        if ( ix(n) /= 0) x(n)=big
      end do
    end do
!
!-----most of the statements in the following section were added 7/8/82

!     check to see if missing data exist or it the
!     data are rectangular (implied missing data).  if so
!     estimate by mean of upper and lower boundaries
!     generated by line of sight method.

    if ( nsim > 3.or.nmiss > 0) go to 7200
    if ( nsim <= 1) go to 3245
!
!  Come here when the data are asymmetric and complete
!  to generate symmetric inicon data.
!
    do i=2,nb
      do j=1,i-1
        xx(i,j)=(xx(i,j)+xx(j,i))/2.0
        xx(j,i)=xx(i,j)
      end do
    end do

    go to 3245
!
!     come here when there is missing data or when the
!     data are rectangulr to flag the missing cells.
!
 7200   n=nn

    do i=nx,nb
      if ( nsim <= 1)i1=i-1
      do j=1,i1
        n=n+1
        if ( ix(n) == 0)go to 135
        xx(i,j)=cut-xx(i,j)
        if ( nsim <= 1)xx(j,i)=xx(i,j)
  135       continue
      end do
    end do

    if ( nsim < 4)go to 130
    ncolp1=ncol+1

    do i=ncol+1,nb
      do j=1,ncol
        xx(j,i)=xx(i,j)
      end do
    end do
!
!  symmetrize the asymmetric matrix
!
  130   if ( debug == 0)go to 3211
    write(lout,3370)
 3370   format(/' data after missing flagged')
    do i=1,nb
      write(lout,6789)(xx(i,j),j=1,nb)
    end do

 3211   if ( nsim <= 1) go to 2010
    na=2
    if ( nsim >= 4)na=ncol+2
    nc=na-1

    do 2011 i=na,nb
    imi=i-1
    do 2011 j=nc,i-1
    if ( xx(i,j) < cut.or.xx(j,i) < cut)go to 2002
 2001   xx(i,j)=(xx(i,j)+xx(j,i))*.5
    xx(j,i)=xx(i,j)
    go to 2011
 2002   if ( xx(i,j) < cut.and.xx(j,i) < cut) go to 2001
    if ( xx(i,j) < cut)xx(i,j)=xx(j,i)
    if ( xx(j,i) < cut)xx(j,i)=xx(i,j)
 2011   continue

 3821   if ( debug == 0)go to 2010
    write(lout,1355)
 1355   format(/' symmetric normalized data')
    do 232 i=1,nb
  232   write(lout,6789)(xx(i,j),j=1,nb)
 6789   format(1x,15f8.2)

!     unless the user has requested, via the noulb parameter,
!     replace missing data esitmates with mean of upper and
!     lower bounds on distances wherever possible.
!     if not possible or if requested use mean calculated above.

 2010   if ( noulb /= 0)go to 1492
    do 150 i=2,nb
    do 150 j=1,i-1
    if ( xx(j,i) > cut)go to 150
    misall=1
    bigdif=0.0
    smlsum=10.0e10
    do 330 ii=1,nb
    if ( ii == i.or.ii == j)go to 330
    ia=i
    ib=j
    ic=ii
    id=ii
    if ( ic >= ia)go to 320
    ia=ii
    ic=i
  320   if ( id >= ib)go to 325
    ib=ii
    id=j
  325   if ( xx(ic,ia) <= cut.or.xx(id,ib) <= cut)go to 330
    misall=0
    sum=    xx(ic,ia)+xx(id,ib)
    dif=abs(xx(ic,ia)-xx(id,ib))
    if ( smlsum > sum)smlsum=sum
    if ( bigdif < dif)bigdif=dif
  330   continue
    if ( misall == 0)xx(j,i)=(smlsum+bigdif)*.5
  150   continue
 1492   continue

    do i=2,nb
      do j=1,i-1
        if ( xx(j,i) < cut)xx(j,i)=cut-xx(j,i)
        xx(i,j)=xx(j,i)
      end do
    end do

    if ( debug == 0)go to 3245
    write(lout,1354)
 1354   format(/' missing data now estimated by mean u/l bounds')
    do 233 i=1,nb
  233   write(lout,6789)(xx(i,j),j=1,nb)
!-----end of added section 7/8/82
!
!  Estimation of additive constant
!
 3245   alph1=0.0
    if ( ndtyp == 1)go to 19
    alph1=xx(2,1)
    do 14 i=2,nb
    imi=i-1
    do 14 j=1,imi
   14   if ( xx(i,j) < alph1) alph1=xx(i,j)
    do 15 i=1,nb
    do 15 j=2,nb
    if ( j == i) go to 15
    j1=j-1
    do 16 k=1,j1
    if ( k == i) go to 16
    tran=xx(i,j)+xx(i,k)-xx(j,k)
    if ( tran < alph1)alph1=tran
   16   continue
   15   continue

   19   wa(1,1,l)=0.0

    do i=2,nb
      wa(i,i,l)=0.0
      imi=i-1
      do j=1,imi
        wa(i,j,l)=(xx(i,j)-alph1)**2
        wa(j,i,l)=wa(i,j,l)
      end do
    end do
!
!-----the following seven lines added 7/8/82
    if ( debug == 0)go to 10
    write(lout,1382)
 1382   format(/' data squared with additive constant estimated', &
        ' (inicon data)')
    do 288 i=1,nb
  288   write(lout,6789)(wa(i,j,l),j=1,nb)
    go to 10
  901   if ( ier == 0)write(lout,9900)
    write(lout,9002)l
9002    format(8x,' all data missing for matrix',i4)
 9900   format(/' alscal fatal error:  computations terminated.')
    write(lout,9001)l,i
    ier=1
   10   nn=nn+nc2

  if ( ier /= 0)return
  rewind ndp
  write(ndp)x
  write(ndp) wa
  rewind ndp

  if ( icnstr == 1 ) then
    nsim = nsim - 2
  end if

  return
end
subroutine step2 ( ix, iy, iz, x, wa, wc, wd, cfr, cfl, w, tr, fk, xx, ws, &
  ds, zz, cv, cw, ndsr, nad, ijkl )
!
!*******************************************************************************
!
!! STEP2 obtains starting configurations and scalings.
!
!
!  Discussion:
!
!    This routine obtains starting configurations and initial optimal
!    scaling to initiate the iterative process
!
!    wa and wc  refer to the same array in the calling routine.
!    this  array  is  passed in as different parameters to permit
!    referencing the arrays using different shape parameters.
!
!  Modified:
!
!    06 September 2002
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  real cfl(nb,1)
  real cfr(nb,1)
  character ( len = 80 ) fmt
  integer ix(nt)
  integer iy(1)
  integer iz(1)
  dimension nad(ns,nb)
  dimension ndsr(ns,nb)
  double precision rnb
  double precision rnbs
  double precision rns
  double precision t
  character ( len = 80 ) title
  double precision trt
  double precision tt
  integer value
  real w(ns,1)

  dimension x(nt),wa(nb,nb,ns),wc(1),wd(1),xx(nb,nb)
  dimension tr(1),fk(ndxp,1)
  dimension ws(nb,1),ds(nb,nb),zz(1,1),cv(1),cw(1)

  double precision cut,stmin
  integer debug,icnstr,noulb
  common /prmblk/cut,stmin,debug,icnstr,noulb
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /block3/ title,fmt
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
  common /inicon/initx,initw,initws,initxc
!
!  compute initial coordinates and weights
!
  istfor=1
  if ( nsim > 3)istfor=2
  rns=1.0 / real ( ns )
  rnb=1.0 / real ( nb )
  rnbs=1.0 / real ( nbs )
  strso=big
  ndim=ndxp-ijkl
!
!     if desired, read initial stimulus coordinates
!
  lll=0
  if ( nsim < 4)go to 3200
  if ( initx >= 2.and.initxc >= 2)go to 3201
  go to 446
 3200 if ( initx < 2)go to 446
 3201 read(in,447)fmt
  447 format(a80)
  k=1
  if ( nsim > 3)k=ncol+1
  do 448 i=k,nb
  448 read(in,fmt)(cfl(i,j),j=1,ndim)
  go to 3026
!
!     otherwise compute  product moment matrices
!
  446 lll=lll+1

  do 24 j=1,nb
  do 24 i=j,nb
   24 xx(i,j)=0.0
  tt=0.0

  do 19 l=1,ns

  trt=0.0
  do 20 i=1,nb
    t=0.0
    do 21 j=1,nb
      t=t+wa(j,i,l)
21  continue
    trt=trt+t
20  tr(i)=t*rnb
  trt=trt*rnbs
  t=0.0
!-----the following two statements added 7/12/82
  if ( debug > 0)write(lout,2144)l
 2144 format(' scalar products for subject',i4)
  do 222 i=2,nb
  do 22 j=1,i-1
    wajil=(wa(j,i,l)-tr(i)-tr(j)+trt)*(-.5)
    wa(j,i,l)=wajil
    wa(i,j,l)=wajil
    t=t+wajil*wajil*2
   22   continue
  waiil=tr(i)-.5*(wa(i,i,l)+trt)
  wa(i,i,l)=waiil
  t=t+waiil*waiil
!-----the following two statements added 7/12/82
  if ( debug > 0)write(lout,2145)(wa(i,j,l),j=1,i)
 2145   format(15f8.2)
  222 continue
  wa11l=tr(1)-.5*(wa(1,1,l)+trt)
  wa(1,1,l)=wa11l
  t=t+wa11l*wa11l
  tt=tt+t
!
!     normalization (within subject)
!
  if ( nsc == 0) go to 19

  if ( debug > 0)write(lout,2143)t
 2143 format(' in step2: t = ',f20.7)
  t=dsqrt(nbs/t)

  do 122 j=1,nb
  do 122 i=j,nb
  wa(i,j,l)=wa(i,j,l)*t
  wa(j,i,l)=wa(i,j,l)
  122 xx(i,j)=xx(i,j)+wa(i,j,l)

19 continue

!     normalization (across subjects)

  if ( nsc /= 0)go to 126
  tt=dsqrt(nbnbns/tt)
  do 123 l=1,ns
  do 123 j=1,nb
  do 123 i=j,nb
  wa(i,j,l)=wa(i,j,l)*tt
  wa(j,i,l)=wa(i,j,l)
  123 xx(i,j)=xx(i,j)+wa(i,j,l)
  126 do 23 j=1,nb
  do 23 i=j,nb
  xx(i,j)=xx(i,j)*rns
   23 xx(j,i)=xx(i,j)
!---nwt and nww are not needed; a later if check is changed
!     nwt=0
!     nww=0
  fmx=0.0
  if ( nwe/2*2 == nwe) go to 25
!
!  compute coordinates and weights for the weighted Euclidean model
!  by using the schonemann-de leeuw method
!
  call init(w,cfl,cfr,nb,ns,ndim,xx,wa,tr,fk,wd,cfl,xx,ws,ndx,ndxp)

  fmx = minval ( w(1:ns,1:ndim) )

  if ( fmx /= 0.0 ) then
    w(1:ns,1:ndim) = w(1:ns,1:ndim) - fmx
  end if

  t=0.0
  do j=1,ndim
    do i=1,ns
      t=t+w(i,j)**2
    end do
  end do

  t=dsqrt((ndim*ns)/t)

  w(1:ns,1:ndim) = t * w(1:ns,1:ndim)

  if ( nwe == 3) go to 4004

  ws(1:nb,1:ndim)=1.0

  go to 26
!
!     compute coordinates for the unweighted Euclidean model
!     by using torgerson's method
!
   25 call cjeig(xx,cfl,cfr,nb,ndim+1,fk,ws ,tr,1,ndxp)
!
!   here ws is being passed to be used as scratch space in cjeig
!
  do 41 j=1,ndim
  tr(j)=sqrt(tr(j))
  do 4080 i=1,nb
  cfl(i,j)=cfl(i,j)*tr(j)
 4080 ws(i,j)=1.0
  do 41 i=1,ns
   41 w(i,j)=1.0
  if ( nwe == 0) go to 26
!
!  Compute coordinates and weights for the asymmetric Euclidean model
!
  cw(1:ndim)=1.0
  go to 4009

 4004 do 4007 j=1,ndim
  t=0.0
  do 4008 i=1,ns
 4008 t=t+w(i,j)
 4007 cw(j)=t*rns
 4009 do 4053 i=1,nb
  ds(i,i)=0.0
  if ( nsim > 1) go to 4053
  do 4052 j=1,i
 4052 ds(j,i)=ds(i,j)
 4053 continue
  call inswm(ds,cfl,cfr,ws,xx,tr,cv,cw,fk,zz,wd,nadct,nb,ndim, &
    (ndim*(ndim+1))/2,ndx,ndxp,ns)
  do 4013 j=1,ndim
  do 4013 i=1,nb
 4013 if ( ws(i,j) < fmx) fmx=ws(i,j)
  if ( fmx >= 0.0)go to 4011
  do 4014 j=1,ndim
  do 4014 i=1,nb
 4014 ws(i,j)=ws(i,j)-fmx
!---nww is not needed
!     nww=1
 4011 t=0.0
  do 4015 j=1,ndim
  do 4015 i=1,nb
 4015 t=t+ws(i,j)**2
  t=dsqrt((ndim*nb)/t)
  do 4016 j=1,ndim
  do 4016 i=1,nb
 4016 ws(i,j)=ws(i,j)*t

!  all paths of the program from above come together to statements 26
!     or 3026

   26 if ( lll == 2)return
!
!     read or print initial configuration and/or initialize
!     initial weights when initial coordinates have been read
!
  if ( initx >= 2)go to 3201
 3026 if ( lll == 2)return
  if ( initxc < 2.or.nsim < 3)go to 3126
  read(in,447)fmt
  do 3030 i=1,ncol
 3030 read(in,fmt)(cfl(i,j),j=1,ndim)
 3126 if ( initw < 2)go to 450
  read(in,447)fmt
  do 449 i=1,ns
  449 read(in,fmt)(w(i,j),j=1,ndim)
  go to 452
  450 if ( initx < 2)go to 452
  do 451 j=1,ndim
  do 451 i=1,ns
  451 w(i,j)=1.0
  452 if ( initws < 2)go to 454
  read(in,447)fmt
  do 453 i=1,nb
  453 read(in,fmt)(ws(i,j),j=1,ndim)
  go to 456
  454 if ( initx < 2)go to 456
  do 455 j=1,ndim
  do 455 i=1,nb
  455 ws(i,j)=1.0
!---initxc added to the following statement 7/6/82
  456 if ( initx > 0.or.initw > 0.or.initws > 0.or.initxc > 0) then
  call page ( lout )
  write(lout,457)
  457   format(' initial configuration'//)
  end if
!---initxc added to the following statement 7/6/82
  if ( initx == 0.and.initxc == 0)go to 461
  write(lout,458)(i,i=1,ndim)
  458 format(' initial stimulus space'/t29,'dimension'/5x,'stimulus',6i12)
  if ( nsim > 3)write(lout,470)
  470 format(/'      column')
  do 459 i=1,nb
  m=i
  if ( nsim > 3.and.i > ncol)m=i-ncol
  ncol1=ncol+1
  if ( nsim > 3.and.i == ncol1)write(lout,471)
  471 format(/' ',8x,'row')
  459 write(lout,460)m,(cfl(i,j),j=1,ndim)
  460 format(i11,5x,6f12.4)
  461 if ( initw == 0)go to 464
  write(lout,462)(i,i=1,ndim)
  462 format(/' initial subject weights'/t29,'dimension'/5x,'subject',6i12)
  do 463 i=1,ns
  463 write(lout,460)i,(w(i,j),j=1,ndim)
  464 if ( initws == 0)go to 467
  write(lout,465)(i,i=1,ndim)
  465 format(/' initial stimulus weights'/t29,'dimension'/5x,'stimulus',6i12)
  do 466 i=1,nb
  466 write(lout,460)i,(ws(i,j),j=1,ndim)
!
!  punch initial configuration if desired
!
  467 if ( nph == 2) then
  write(nplt,3001)
 3001   format('(6f13.9)')
  do 3002 i=1,nb
 3002   write(nplt,3003)(cfl(i,j),j=1,ndim)
 3003   format(6f13.9)
  if ( nwe/2*2 == nwe)go to 3005
  write(nplt,3001)
  do 3004 i=1,ns
 3004   write(nplt,3003)(w(i,j),j=1,ndim)
 3005   if ( nwe >= 2) then
    write(nplt,3001)
    do 3006 i=1,nb
 3006     write(nplt,3003)(ws(i,j),j=1,ndim)
  end if
  end if
!
!  Prepare qualitative data for optimal scaling.
!
  call page ( lout )
  write(lout,320) title,ndim
  320 format(a80//' iteration history for the',i2,'  dimensional solution')
!-pc----------------------
  call clear
  write(*,321) istfor
  write(lout,321) istfor
  321 format(' sstress (in squared distances) formula',i2,' is used.'/// &
  6x,'iteration',6x,'sstress',6x,'improvement'/)

  if ( ijkl > 1)go to 30
  iz(1)=1
  n2=0
  if ( ndtyp-3) 30,31,144
   31 if ( nwc == 1) go to 431
  if ( nwc == 2) go to 2012
!
!     unconditional ordinal data
!
  do k=1,nt
    ix(k)=k
    wc(k)=x(k)
  end do

  call shel9(wc,ix,nt)

  value = iz(1)
  call bloc2 ( wc, iy, value, nt )

  go to 30
!
!     matrix conditional ordinal data
!
  431 continue

  do l=1,ns

    do k=1,nc2
      m=n2+k
      ix(m)=k
      wc(m)=x(m)
    end do

    n=n2+1
    call shel9(wc(n),ix(n),nc2)

    if ( l /= 1 ) then
      iz(l)=iz(l-1)
    end if

    call bloc2(wc(n),iy,iz(l),nc2)
    n2=n2+nc2

  end do

  go to 30
!
!     row conditional ordinal data
!
 2012 continue

  n1=0

  do l = 1, ns

    do i = 1, nb

      do j = 1, nb
        m = n2+j
        ix(m) = j
        wc(m) = x(m)
      end do

      n=n2+1
      call shel9(wc(n),ix(n),nb)
      n1=n1+1

      if ( n1 /= 1 ) then
        iz(n1)=iz(n1-1)
      end if

      call bloc2(wc(n),iy,iz(n1),nb)
      n2=n2+nb

    end do

  end do

  go to 30

  144 continue

  if ( nwc == 1) go to 573
  if ( nwc == 2) go to 2017
!
!     unconditional nominal data
!
  ndct=0
  do 145 l=1,nt
  if ( ix(l) == 1) go to 1145
  if ( ndct == 0) go to 1146
  do 146 k=1,ndct
  if ( wc(k) == x(l)) go to 147
  146 continue
 1146 ndct=ndct+1
  wc(ndct)=x(l)
  ix(l)=ndct
  go to 145
  147 ix(l)=k
  go to 145
 1145 ix(l)=-1
  145 continue
  if ( nab == nt) go to 30
  ndct=ndct+1
  do 2310 l=1,nt
 2310 if ( ix(l) == -1) ix(l)=ndct
  go to 30

!     matrix conditional nominal data

  573 do 574 l=1,ns
  ndsbb=0
  n1=n2+1
  n2=n2+nc2
  do 575 j=n1,n2
  if ( ix(j) == 1) go to 1575
  if ( ndsbb == 0) go to 1577
  do 576 k=1,ndsbb
  if ( wc(k) == x(j)) go to 577
  576 continue
 1577 ndsbb=ndsbb+1
  wc(ndsbb)=x(j)
  ix(j)=ndsbb
  go to 575
  577 ix(j)=k
  go to 575
 1575 ix(j)=-1
  575 continue
  if ( nad(l,1) == nc2) go to 574
  ndsbb=ndsbb+1
  do 2311 j=n1,n2
 2311 if ( ix(j) == -1) ix(j)=ndsbb
  574 ndsr(l,1)=ndsbb
  go to 30
!
!     row conditional nominal data
!
 2017 do 2018 l=1,ns
  do 2019 i=1,nb
  ndsrr=0
  n1=n2+1
  n2=n2+nb
  do 2020 j=n1,n2
  if ( ix(j) == 1) go to 2021
  if ( ndsrr == 0) go to 2022
  do 2023 k=1,ndsrr
  if ( wc(k) == x(j)) go to 2024
 2023 continue
 2022 ndsrr=ndsrr+1
  wc(ndsrr)=x(j)
  ix(j)=ndsrr
  go to 2020
 2024 ix(j)=k
  go to 2020
 2021 ix(j)=-1
 2020 continue
  if ( nad(l,i) == nb) go to 2019
  ndsrr=ndsrr+1
  do 2025 j=n1,n2
 2025 if ( ix(j) == -1) ix(j)=ndsrr
 2019 ndsr(l,i)=ndsrr
 2018 continue
   30 continue
!
!     perform initial optimal scaling
!
!---the following two statements are replaced by the next statement
!     if ( nwe == 0.and.ndtyp /= 1) go to 8888
!     if ( nwt == 0.and.nww == 0)return
  if ( nwe == 0.and.ndtyp == 1)return
 8888 call distp(w,cfl,x,wc,wd,ix,iy,iz,xx,ndsr,ws,nad)
  if ( istfor == 2)return
  if ( nwe == 0.and.strss < 0.5)return
  if ( initx > 1.or.initw > 1.or.initws > 1.or.initxc > 1)return
  write(lout,8889)strss
  write(*,8889)strss
 8889 format(10x,'0',2x,f15.5)
  n=0
  do 442 l=1,ns
  if ( nsim > 1) go to 2026
  do 444 i=1,nb
  444 wa(i,i,l)=0.0
  do 443 i=2,nb
  i1=i-1
  do 443 j=1,i1
  n=n+1
  wa(i,j,l)=x(n)
  wa(j,i,l)=x(n)
  443 continue
  go to 442
 2026 do 2027 i=1,nb
  do 2027 j=1,nb
  n=n+1
 2027 wa(i,j,l)=x(n)
  do 2028 i=1,nb
  do 2028 j=1,i
  wa(i,j,l)=(wa(i,j,l)+wa(j,i,l))*.5
  wa(j,i,l)=wa(i,j,l)
 2028 continue
  442 continue
  go to 446

end
subroutine step3 ( ix, iy, iz, x, g, wa, wb, wc, u11, u12, u22, ub1, ub2, xn, &
  xx, cfl, w, tr, fk, ws, ndsr, nad, ijkl, iret )
!
!*******************************************************************************
!
!! STEP3 is the major computational routine.
!
!
!  Discussion:
!
!    this is the major computational routine. it controls the flow
!    of the iterative minimization process, performs regressions
!    to calculate weights, and outputs the results
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!    wa,wb and wc refer to the same array in the calling routine
!    these arrays are passed in as different parameters to permit
!    referencing the arrays using different numbers of dimensions
!
  double precision t,tt,a,b,c
  integer ix(nt),iy(1),iz(1)
  dimension x(nt),g(nbs,1),wa(nb,nb,ns),wb(nbs,ns),wc(1)
  dimension u11(ndx,1),u12(ndx,1),u22(ndx,1),ub1(1),ub2(1),xn(1)
  dimension xx(nb,nb),cfl(nb,1),w(ns,1),tr(1),fk(ndxp,1)
  dimension ws(nb,1),ndsr(ns,nb),nad(ns,nb)
  double precision cut,stmin
  integer debug,icnstr,noulb
  common /prmblk/cut,stmin,debug,icnstr,noulb
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
  common /inicon/initx,initw,initws,initxc
!
!     iteration begins here
!
  iret=0
  do 2 ll=1,maxit
!
!  Compute squared distances,disparities and stress
!
  call distp(w,cfl,x,wc,g,ix,iy,iz,xx,ndsr,ws,nad)

  dif=strso-strss

  if ( ll == 1) then
    write(lout,321)ll,strss
    write(*,321)ll,strss
  else
    write(lout,321)ll,strss,dif
    write(*,321)ll,strss,dif
  end if

  321 format(i11,2x,2f15.5)
!     test for exit
  if ( dif < 0.0) go to 145
!----- following statement added 7/7/82 (sas 10/23/81)
  if ( strss < stmin) go to 141
  if ( dif < epsi) go to 143
  if ( ll >= maxit) go to 3
  strso=strss
!
!     set-up for the next iteration
!
  a = 1.0 / (1.0-strss*strss)
  if ( nsim > 1) go to 2029
!   nsim <= 1 implies nc2=nb*(nb-1)/2 therefore
  n=0
  do 43 l=1,ns
  wa(1,1,l)=0.0
  do 44 i=2,nb
  i1=i-1
  wa(i,i,l)=0.0
  do 44 j=1,i1
  n=n+1
  wa(i,j,l)=x(n)*a
   44 wa(j,i,l)=wa(i,j,l)
   43 continue
  go to 2030
 2029 do 2031 l=1,nbnbns
 2031 wc(l)=x(l)*a
 2030 continue
  if ( nwe == 0) go to 4017
  l=0
  do 47 i=1,nb
  do 47 j=1,nb
  l=l+1
  do 47 k=1,ndim
   47 g(l ,k)=(cfl(i,k)-cfl(j,k))**2
  if ( nwe == 2) go to 52
  if ( initw == 3)go to 4533

!     (1) estimate subject weights

  do 48 i=1,ndim
  do 48 j=1,i
  t=0.0
  l=0
  do 49 kk=1,nb
  tt=ws(kk,i)*ws(kk,j)
  do 49 k=1,nb
  l=l+1
   49 t=t+g(l,i)*g(l,j)*tt
  fk(i,j)=t
   48 fk(j,i)=t
  call minv(fk,ndim,ndxp)
  do 51 j=1,ndim
  do 51 i=1,ns
  t=0.0
  do 50 k=1,ndim
  l=0
  do 50 m=1,nb
  tt=ws(m,k)*fk(j,k)
  do 50 n=1,nb
  l=l+1
   50 t=t+wb(l,i)*g(l,k)*tt
   51 w(i,j)=t
!
!     check for negative subject weights
!
  if ( nnc == 1) go to 4533
  do 471 k=1,ns
  do 476 j=1,ndim
  476 tr(j)=w(k,j)
  do 472 j=1,ndim
  if ( tr(j) < 0.0) go to 473
  472 continue
  go to 471
  473 tr(j)=0.0
  if ( ndim == 1) go to 4070
  j1=j+1
  if ( j == ndim) j1=1
  do 482 lx=1,20
  do 478 i=j1,ndim
  a=0.0
  b=0.0
  do 477 i1=1,nbs
  kk=(i1-1)/nb+1
  c=wb(i1,k)
  do 474 i2=1,ndim
  474 if ( i2 /= i)c=c-g(i1,i2)*ws(kk,i2)*tr(i2)
  a=a+c*g(i1,i)*ws(kk,i)
  477 b=b+g(i1,i)**2*ws(kk,i)**2
  tr(i)=a/b
  478 if ( tr(i) < 0.0) tr(i)=0.0
  do 479 j=1,ndim
  if ( abs(tr(j)-w(k,j)) > 0.0005) go to 483
  479 continue
  go to 4070
  483 do 481 j=1,ndim
  481 w(k,j)=tr(j)
  j1=1
  482 continue
  go to 471
 4070 do 484 j=1,ndim
  484 w(k,j)=tr(j)
  471 continue
  do 355 j=1,ndim
  do 356 i=1,ns
  if ( w(i,j) /= 0.0) go to 355
  356 continue
  go to 4
  355 continue
 4533 if ( nwe == 1) go to 4017
   52 if ( initws == 3)go to 4017

!     (2) estimate stimulus weights

  do 4022 k=1,nb
  do 4018 i=1,ndim
  do 4018 j=1,i
  t=0.0
  n=(k-1)*nb
  do 4019 l=1,nb
  n=n+1
  tt=g(n,i)*g(n,j)
  do 4019 kk=1,ns
 4019 t=t+tt*w(kk,i)*w(kk,j)
  fk(j,i)=t
 4018 fk(i,j)=t
  call minv(fk,ndim,ndxp)
  do 4021 j=1,ndim
  t=0.0
  do 4020 i=1,ndim
  kk=(k-1)*nb
  do 4020 n=1,nb
  kk=kk+1
  tt=g(kk,i)*fk(i,j)
  do 4020 l=1,ns
 4020 t=t+tt*wb(kk,l)*w(l,i)
 4021 ws(k,j)=t

!     check for negative stimulus wieghts

  if ( nnc == 1) go to 4022
  do 4026 j=1,ndim
 4026 tr(j)=ws(k,j)
  do 4023 j=1,ndim
  if ( tr(j) < 0.0) go to 4024
 4023 continue
  go to 4022
 4024 tr(j)=0.0
  if ( ndim == 1) go to 4071
  j1=j+1
  if ( j == ndim) j1=1
  do 4025 lx=1,20
  do 4027 i=j1,ndim
  a=0.0
  b=0.0
  do 5029 i1=1,nb
  ii1=(k-1)*nb+i1
  do 5029 i11=1,ns
  c=wb(ii1,i11)
  do 4028 i2=1,ndim
 4028 if ( i2 /= i)c=c-g(ii1,i2)*tr(i2)*w(i11,i2)
  a=a+c*g(ii1,i)*w(i11,i)
 5029 b=b+g(ii1,i)**2*w(i11,i)**2
  tr(i)=a/b
 4027 if ( tr(i) < 0.0) tr(i)=0.0
  do 4030 j=1,ndim
  if ( abs(tr(j)-ws(k,j)) > 0.0005) go to 4031
 4030 continue
  go to 4071
 4031 do 4033 j=1,ndim
 4033 ws(k,j)=tr(j)
  j1=1
 4025 continue
  go to 4022
 4071 do 4032 j=1,ndim
 4032 ws(k,j)=tr(j)
 4022 continue
  do 4034 j=1,ndim
  do 4035 i=1,nb
  if ( ws(i,j) /= 0.0) go to 4034
 4035 continue
  go to 4
 4034 continue

!     (3) estimate configuration

 4017 if ( initx == 3.and.nsim < 4)go to 2
  if ( initx == 3.and.initxc == 3.and.nsim >= 4)go to 2
  call inner(cfl,w,x,wb,u11,u12,u22,ub1,ub2,xn,nb,ndim,ns,ndx,nbs,ws)
    2 continue

    3 write(lout,9901)
  write(*,9901)
 9901 format(/' alscal message:  iterations stopped because ')
  write(lout,299)maxit
  write(*,299)maxit
  299 format(t19,'this is iteration',i4//)
  strss=strso
  go to 9999
!----- the following four statements added 7/7/82 (sas 10/23/81)
  141 write(lout,9901)
  write(*,9901)
  write(lout,142)stmin
  write(*,142)stmin
  142 format(t19,'s-stress less than',f9.6//)
  go to 9999
  143 write(lout,9901)
  write(*,9901)
  write(lout,144)epsi
  write(*,144)epsi
  144 format(t19,'sstress improvement less than',f9.6//)
  go to 9999

  145 continue
  write ( *, * ) ' '
  write ( *, * ) 'ALSCAL - Warning!'
  write ( *, * ) '  Iteration terminated, negative stress improvement.'
  write ( *, * ) '  Data may be ill-conditioned,'
  write ( *, * ) '  or there may be a program error.'
  go to 9999

    4 write(lout,9900)
  write(*,9900)
 9900 format(/' alscal warning:    computations interrupted.')
  write(lout,298)
  write(*,298)
  298 format(8x,'a dimension has only zero weights. solution skipped')
  if ( ijkl == nd) go to 4099
  if ( npt == 1.and.nt == nab) go to 2899
  rewind ndp
  read(ndp) x
 2899 read(ndp) wa
 4099 rewind ndp
  iret=1
  return

 9999 call hitr

  return
end
subroutine step3a ( ix, iy, iz, x, wc, wd, xx, cfl, w, ws, ndsr, nad, phisub, &
  phisti )
!
!*******************************************************************************
!
!! STEP3A prints stress and correlations in distances.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
!     subroutine to print stress and correlations in distances
!        normalizes weights according to correlation
!
  integer ix(nt)
  integer iy(1)
  integer iz(1)
  dimension x(nt),wc(1),wd(1),xx(nb,nb),cfl(nb,1)
  dimension w(ns,1),ws(nb,1),phisub(ns),phisti(nb)
  dimension ndsr(ns,nb),nad(ns,nb)
  character ( len = 120 ) line
!
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata

  line=' '
  istfor=1
  if ( nsim > 3)istfor=2
!
!  stress in distances
!
  do i=1, nt

    if ( x(i) >= 0.0 ) then
      x(i)=sqrt(x(i))
    else
      x(i)=-sqrt(-x(i))
    end if

    if ( wc(i) >= 0.0 ) then
      wc(i)=sqrt(wc(i))
    else
      wc(i)=-sqrt(-wc(i))
    end if

  end do

  write(lout,8451)
 8451 format(////' stress and squared correlation (rsq) in distances'/ &
     /' rsq values are the proportion of variance of the scaled data ', &
     '(disparities) in the partition'/' (row, matrix, or entire data)', &
     ' which is accounted for by their corresponding distances.')
  write(lout,7483)istfor

  call clear
  write (*,8452)
 8452 format(/' stress and squared correlation (rsq) in distances'/)

 7483 format(/' kruskal''s stress formula',i2,' is used.'/)
  if ( nwc-1)8460,8470,8480
!
!  unconditional stress output
!
 8460 call mstrs(nt,x,wc,strss1,strss2)
  phi=1.0-strss2
  strss1=sqrt(strss1)
  write(lout,8461)strss1,phi
  write(*,8461)strss1,phi
 8461 format(' stress = '    ,f5.3,'      rsq = ',f5.3)
  go to 8490
!
!  matrix conditional stress output
!
 8470 stres1=0.0
  stres2=0.0
  if ( ns /= 1)write(lout,8473)
 8473 format(4('   matrix     stress      rsq '))
  n=1
  do 8479 l=1,ns
  call mstrs(nc2,x(n),wc(n),strss1,strss2)
  stres1=stres1+strss1
  stres2=stres2+strss2
  phisub(l)=1.0-strss2
  strss1=sqrt(strss1)
  lmod4=mod(l-1,4)
  lstart=5+lmod4*30
  lend=lstart+24
  if ( ns /= 1)write(line(lstart:lend),8471)l,strss1,phisub(l)
 8471 format(i3,f13.3,f9.3)
  if ( lmod4 == 3)write(lout,8475)line
 8475 format(a)
  n=n+nc2
 8479 continue
  if ( lmod4 /= 3)write(lout,8475)line(1:29+lmod4*30)
  line=' '
  phi=1.0-stres2/ns
  stres1=sqrt(stres1/ns)
  if ( ns /= 1)write(lout,8472)stres1,phi
  if ( ns /= 1)write(lout,9472)stres1,phi
 8472 format(/'   overall',f10.3, f9.3/)
  if ( ns == 1)write(lout,8461)stres1,phi
  if ( ns == 1)write(*,8461)stres1,phi
  go to 8490
!
!  row conditional stress output
!
 8480 stres1=0.0
  stres2=0.0
  do 8485 j=1,nb
 8485 phisti(j)=0.0
  n1=1

  nbornr=nb
  if ( nsim > 3)nbornr=nb-ncol

  do 8489 i=1,ns
  write(lout,8481) i
 8481 format(' matrix',i4/)
  if ( nsim > 3)write(lout,2000)
 2000 format(' row stimuli')
  write(lout,2002)
 2002 format(4('  stimulus    stress      rsq '))
  str1=0.0
  str2=0.0
  n=n1
  do 8488 j=1,nb

!----- the following statement was moved to correct and error in
!----- calculating average stress and rsq indices for rectangular
!----- row conditional data. 7/6/82 (sas change 12/1/80)
!----- the error was that the perfect fit indices for columns
!----- were being included in the summation process.

  if ( nsim > 3.and.j <= ncol)go to 2001
  call mstrs(nb,x(n), wc(n), strss1,strss2)
  str1=str1+strss1
  str2=str2+strss2
  phirow=1.0-strss2
  phisti(j)=phisti(j)+phirow
  strss1=sqrt(strss1)
!     if ( nsim > 3.and.j <= ncol)go to 2001
  k=j
  if ( nsim > 3)k=j-ncol
  lmod4=mod(k-1,4)
  lstart=5+lmod4*30
  lend=lstart+24
  write(line(lstart:lend),8471)k,strss1,phirow
  if ( lmod4 == 3)write(lout,8475)line
 2001 n=n+nb
 8488 continue
  if ( lmod4 /= 3)write(lout,8475)line(1:29+lmod4*30)
  line=' '
  stres1=stres1+str1
  stres2=stres2+str2

!----- the following two statements were changed to compute
!----- average stress and rsq indices correctly for rectangular
!----- row conditional data.  the error was that the division
!----- was by the total number of simuli, not the number of rows.
!----- 7/6/82 (sas change 12/1/80)
!     phisub(i)=1.0-str2/nb
!     str1=sqrt(str1/nb)

  phisub(i)=1.0-str2/nbornr
  str1=sqrt(str1/nbornr)
  write(lout,8472)str1,phisub(i)
  n1=n1+nc2
 8489 continue
  if (ns == 1)go to 8490
  write(lout,8491)
 8491 format(/' averaged over matrices')
  if ( nsim < 4)write(lout,8492)
  if ( nsim > 3)write(lout,8493)
 8492 format(/'  stimulus',6x,'rsq')
 8493 format(/'  row stimulus  rsq')
  do 8486 i=1,nb
  phisti(i)=phisti(i)/ns
  k=i
  if ( nsim > 3.and.k <= ncol)go to 8486
  if ( nsim > 3)k=i-ncol
  write(lout,8487)k,phisti(i)
 8486 continue
 8487 format(t5,i3,f13.3,f9.3)
  na=nb
  if ( nsim > 3)na=nb-ncol
  phi=1.0-stres2/(na*ns)
  stres1=sqrt(stres1/(na*ns))
  write(lout,9472)stres1,phi
  write(*,9472)stres1,phi
 9472 format(/' overall stress =',f6.3,' and rsq =',f6.3)
!
!     normalize derived solution
!
 8490 call normx(cfl,w,nb,ns,ndim,nwe)

  iconfl = 0
  if ( nwc == 0 )iconfl=1

  if ( nwe == 1 .or. nwe == 3 ) then
    call normw(w,ns,ndim,phi,phisub,iconfl)
  end if

  iconfl=1
  if ( nwc == 2)iconfl=0
  if (nwe >= 2) call normw(ws,nb,ndim,phi,phisti,iconfl)
  nad(1,1)=-nad(1,1)
  call distp(w,cfl,x,wc,wd,ix,iy,iz,xx,ndsr,ws,nad)

  call hitr

  return
end
subroutine step4 ( x, wc, disp, cfl, w, tr, ws, dist, dummy, pijp, xk, row, &
  ijkl )
!
!*******************************************************************************
!
!! STEP4 outputs the results.
!
!
!  Author:
!
!    Copyright, 1977, 
!    Forrest W. Young, Yoshio Takane and Rostyslaw Lewyckyj
!
  integer nb
  integer ns
  integer nt
!
  character alfa(52)
  character ( len = 52 ) alfa1
  character alfal(26)
  character alfas(26)
  real cfl(nb,1)
  real disp(nb,nb)
  character ( len = 80 ) fmt
  character ( len = 80 ) title
  real tr(1)
  real w(ns,1)
  real wc(nbnbns)
  real wl(6)
  real ws(nb,1)
  real x(nt)
!
  equivalence(alfa,alfa1),(alfa1(1:26),alfal),(alfa1(27:52),alfas)
!-----following three statements added for gemscal 8/9/82
  dimension dist(nb,nb),dummy(nb,nb),pijp(nb,nb),xk(nb,6),row(nb)
  common /pdcom/iflpds,ndir
  common /block1/nc,nd,big,nc2,ndt,nnc,nph,npt,nsc, &
         epsi,ndim,ndx,ndxs,ndxp,maxit,nadct,ndct,strso, &
         strss,strss2,nb,ns,ndtyp,nps,nwc,ndeg,nt,nbs,nbnbns
  common /block2/ncst,nsim,nwe,ndmx,nab,ncol
  common /block3/ title,fmt
  common /ionums/in,nplt,lout,ndp,ndq,ndr,ndpp,indata
!-pc--------------------------------------------------------------------
  data alfa1/'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'/
!-----------------------------------------------------------------------

!-----following four statements added 8/9/82
  if ( ndir == 0)ndir=ndim
!-----following statement added 19/10/82 for gemscal
  iplot=npt
  do 80 i=1,nt
  if ( x(i) < 0.0) go to 4081
  x(i)=sqrt(x(i))
  go to 4082
 4081 x(i)=-sqrt(-x(i))
 4082 if ( wc(i) < 0.0) go to 4083
  wc(i)=sqrt(wc(i))
  go to 80
 4083 wc(i)=-sqrt(-wc(i))
   80 continue
!
!  Print the results.
!
  if ( ndim > 1 ) then
    call arnge ( w, cfl, ws, nb, ns, ndim, tr )
  end if

  write(lout,272) ndim
  272 format(////' configuration derived in',i3,'  dimensions')
  write(lout,228) (j,j=1,ndim)
  228 format(///' stimulus coordinates'/27x,'dimension'/ &
     5x,'stimulus',2x,'plot',i6,5i12)
  write(lout,2683)
 2683 format(6x,'number',2x,'symbol')
  write(lout,54)

  if ( nsim > 3) then
  write(lout,2684)
 2684   format(6x,'column')
    do 6401 i=1,ncol
      kk=mod(i-1,26)+1
 6401     write(lout,521) i,alfas(kk),(cfl(i,j),j=1,ndim)
  521     format(i10,7x,a1,f10.4,5f12.4)
  write(lout,2685)
 2685   format(/,8x,'row')
    do 6402 i=ncol+1,nb
      m=i-ncol
      kk=mod(m-1,26)+1
 6402     write(lout,521) m,alfal(kk),(cfl(i,j),j=1,ndim)
  else
  do 64 i=1,nb
    kk=mod(i-1,52)+1
   64   write(lout,521) i,alfa(kk),(cfl(i,j),j=1,ndim)
  end if

  if ( (nwe/2)*2 == nwe) go to 62
!-----most of the statements in the following section added 8/9/82
  if ( nwc /= 0)write(lout,4200)
 4200 format(//'subject weights measure the importance of each dimension', &
      ' to each subject.  squared weights sum to rsq.'// &
      'a subject with weights proportional to the average ' &
      'weights has a weirdness of zero, the minimum value.'/ &
      'a subject with one large weight and many low weights has ', &
      'a weirdness near one.'/ &
      'a subject with exactly one positive weight has a weirdness ' &
      'of one, the maximum value for nonnegative weights.')
  write(lout,224) (j,j=1,ndim)
  224 format(//' subject weights'/39x,'dimension'/ &
      'subject   plot    weird-',i6,5i12)
  write(lout,1050)
 1050 format(6x,'number  symbol    ness')
  write(lout,54)
   54 format(' ')
!
!  Calculate weirdness.
!
  do j=1,ndim
    wl(j)=0.0
    do i=1,ns
      wl(j)=wl(j)+w(i,j)
    end do
  end do

  dimnum = ndim
  cosmax = 1.0 / sqrt(dimnum)

  do i = 1, ns

    top=0.0
    bot=0.0
    do j=1,ndim
      val=w(i,j)/wl(j)
      top=top+val
      bot=bot+val*val
    end do

    weird=acos(cosmax*top/sqrt(bot))/acos(cosmax)
    kk=mod(i-1,52)+1
    write(lout,521)i,alfa(kk),weird,(w(i,j),j=1,ndim)

  end do

  do j=1,ndim
    wl(j)=0.0
    do i=1,ns
      wl(j)=wl(j)+w(i,j)**2
    end do
    wl(j)=wl(j)/ns
  end do

  write(lout,1004)(wl(j),j=1,ndim)
 1004 format(/'     average (rms)',10x,6f12.4)

   62 if ( nwe < 2) go to 4040
  write(lout,352) (j,j=1,ndim)
  352 format(//' stimulus weights'/27x,'dimension'/5x,'stimulus', &
      2x,'plot',i6,5i12)
  write(lout,2683)
  write(lout,54)
  do 4042 i=1,nb
  kk=mod(i-1,52)+1
 4042 write(lout,521)i,alfa(kk),(ws(i,j),j=1,ndim)
 4040 continue
!
!     plot results
!
 3321 if ( npt == 0) go to 4043
  xhi=2.5
  xlo=-2.5
  yhi=2.5
  ylo=-2.5
  if ( ndim == 1) go to 66

  do i=2,ndim
    i1=i-1
    do j=1,i-1
      call page ( lout )
      write(lout,229) title,j,i
  229     format(a80/6x,'derived stimulus configuration',/6x, &
          'dimension',i3,' (horizontal)  vs  dimension',i3,' (vertical)')
      call plotr(cfl(1,j),cfl(1,i),xhi,yhi,xlo,ylo,nb,lout,2,nb)
    end do
  end do

  go to 67
   66 call page ( lout )
  write(lout,230) title
  230 format(a80,/6x,'derived stimulus configuration',/6x, &
      'one dimensional plot')
  call plotr(cfl(1,1),cfl(1,1),xhi,yhi,xlo,ylo,nb,lout,2,nb)
!-----code relating to fmw,fmz,fmu and fmv changed 8/10/82
   67 fmw=0.9
  fmz=0.0
  if ( nnc == 1)fmz=-fmw
  fmu=1.0
  fmv=0.0
  if ( nnc == 1)fmv=-fmu
  iax=2
  if ( nnc /= 1)iax= -2
  if ( nwe/2*2 == nwe) go to 70

  do i = 2, ndim
    do j=1,i-1
      call page ( lout )
      write(lout,231) title,j,i
  231     format(a80/6x,'derived subject weights'/6x, &
          'dimension',i3,' (horizontal)  vs  dimension',i3,' (vertical)')
      call plotr(w(1,j),w(1,i),fmu,fmw,fmv,fmz,ns,lout,iax,ns)
    end do
  end do

   70 if ( nwe < 2) go to 4043
!-----code relating to fmu,...,iax replaced by code above 8/10/82
  do 4046 i=2,ndim
  i1=i-1
  do 4046 j=1,i1
  call page ( lout )
  write(lout,4047)title,j,i
 4047 format(a80,/6x,'derived stimulus weights'/6x, &
      'dimension',i3,' (horizontal)  vs  dimension',i3,' (vertical)')
 4046 call plotr(ws(1,j),ws(1,i),fmu,fmw,fmv,fmz,nb,lout,iax,nb)
 4043 continue
!
!  Fit generalized Euclidean model weights
!
  if ( ndt /= 1.and.iflpds /= 1)go to 5099
  knt=1
  do 5006 i=1,ns
  if ( ndt /= 1)go to 5010
  call page ( lout )
  write(lout,5001)title,i
 5001 format(a80,/6x,'optimally scaled data (disparities)     subject',i4)
  if ( nsim < 2)call outs(x(knt),nb,1,1,disp)
  if ( nsim > 1)call outa(x(knt),nb,-1)
 5010 if ( iflpds /= 1)go to 5005
  call page ( lout )
  write(lout,5002)title,i
 5002 format(a80/1x,'general Euclidean model:    subject',i4)
  if ( nsim > 1)go to 5004
!
!  Symmetric data
!
  disp(1,1)=0.0

  do j = 2, nb
    disp(j,j) = 0.0
    do k = 1, j-1
      disp(k,j) = x(knt)
      disp(j,k) = x(knt)
      knt = knt + 1
    end do
  end do

  call pdmain(cfl,disp,nb,ndim,ndir,i,nsim,iplot,lout,dist,dummy,pijp,xk,row)

  go to 5006
!
! asymmetric and rectangular data
!
 5004 continue
  call pdmain(cfl,x(knt),nb,ndim,ndir,i,nsim,iplot,lout, &
    dist,dummy,pijp,xk,row)
 5005 if ( nsim < 2)knt=knt+nb*(nb-1)/2
  if ( nsim > 1)knt=knt+nb*nb

 5006 continue
!-----end of gemscal section added 8/10/82
 5099 if ( npt == 0)go to 5265
  call page ( lout )
  write(lout,235) title
  235 format(a80/6x,'scattergram (plot of linear fit)'/ &
     6x,'distances (vertical) vs disparities (horizontal)')
  rewind ndq
  write(ndq) x
  nn=nt
  if ( nsim < 4)go to 4065
  l=0
  ll=0
  do 4064 k=1,ns
  do 4064 i=1,nb
  do 4064 j=1,nb
  l=l+1
  if ( i <= ncol)go to 4064
  if ( j > ncol)go to 4064
  ll=ll+1
  x(ll)=x(l)
  wc(ll)=wc(l)
 4064 continue
  nn=ll
 4065 call plotr(x,wc,1.0,1.0,1.0,1.0,nn,lout,-1,nt)
  rewind ndq
  rewind ndpp
  i1=nb
  nx=2
  if ( nsim > 1)nx=1

  l=0
  do k=1,ns
    do i=nx,nb
      if ( nsim <= 1)i1=i-1
      i2=l+1
      l=l+i1
      read(ndpp)(x(ii),ii=i2,l)
    end do
  end do

  if ( ndeg == 1.and.ndtyp <= 2)go to 5265
  if ( nsim < 4)go to  4075
  l=0
  ll=0
  do 4074 k=1,ns
  do 4074 i=1,nb
  do 4074 j=1,nb
  l=l+1
  if ( i <= ncol)go to 4074
  if ( j > ncol)go to 4074
  ll=ll+1
  x(ll)=x(l)
 4074 continue
 4075 if ( nwc == 2)go to 4076
  if ( nwc == 1.and.ns > 1)go to 4076
  call page ( lout )
  write(lout,233) title
  233 format(a80/6x,'plot of nonlinear fit'/ &
     6x,'distances (vertical) vs observations (horizontal)')
  call plotr(x,wc,1.0,1.0,1.0,1.0,nn, lout,-1,nt)
 4076 read(ndq) (wc(i),i=1,nt)
  rewind ndq
  if ( nsim < 4)go to 4085
  l=0
  ll=0
  do 4084 k=1,ns
  do 4084 i=1,nb
  do 4084 j=1,nb
  l=l+1
  if ( i <= ncol)go to 4084
  if ( j > ncol)go to 4084
  ll=ll+1
  wc(ll)=wc(l)
 4084 continue
!-----plotting of transformations updated 8/10/82
 4085 if ( nwc == 2)go to 5100
  if ( nwc == 1.and.ns > 1)go to 5100
!
!  Make one plot when there is one nonlinear transformation
!
  call page ( lout )
  write(lout,234) title
  234 format(a80/6x,'plot of transformation'/ &
      6x,'disparities (vertical) vs observations (horizontal)')
  call plotr(x,wc,1.0,1.0,1.0,1.0,nn, lout,-1,nt)
 5100 if ( iplot <= 1.or.nwc == 0)go to 5265
  if ( nwc == 1.and.ns == 1)  go to 5265
  loc=1
  if ( nwc == 2)go to 5110

! make several plots when fullplot requested
! and there are several transformations

  do 5102 k=1,ns
  call page ( lout )
  write(lout,5101)title,k
 5101 format(a80/6x,'plot of transformation     subject',i4/ &
      6x,'disparities (vertical) vs observations (horizontal)')
  call plotr(x(loc),wc(loc),1.0,1.0,1.0,1.0,nc2,lout,-1,nc2)
 5102 loc=loc+nc2
  go to 5265

! do fullplot when row conditional

 5110 do 5112 k=1,ns
  do 5112 i=1,nb
  call page ( lout )
  write(lout,5111)title,k,i
 5111 format(a80/6x,'plot of transformation     subject',i4, &
  ', row',i4/6x,'disparitites (vertical) vs observations (horizontal)')
  call plotr(x(loc),wc(loc),1.0,1.0,1.0,1.0,nb,lout,-1,nb)
 5112 loc=loc+nb

!     punch results

 5265 if ( nph >= 1) then
  write(nplt,252)
  do 71 i=1,nb
   71   write(nplt,251) (cfl(i,j),j=1,ndim)
  251   format(6f13.9)
  252   format('(6f13.9)')
  if ( nwe/2*2 /= nwe) then
    write(nplt,252)
    do i=1,ns
      write(nplt,251) (w(i,j),j=1,ndim)
    end do
  end if

  if ( nwe >= 2) then
    write(nplt,252)
    do i=1,nb
      write(nplt,251) (ws(i,j),j=1,ndim)
    end do
  end if

  end if

!-pc----------------------------------------------------
! write file "cplot.dat" for rotation of configuration
  if (nph == 1) call wplot (cfl,w,ws,nwe,ndim,nb,ns)
! ------------------------------------------------------

  if ( ijkl == nd) go to 4099
  rewind ndp
  read(ndp) x
  read(ndp) wc
 4099 rewind ndp
!-----flattened subject weight section added 8/10/82
  if ( nwe/2*2 == nwe)go to 9999
  do 4102 i=1,ns
  sum=0.0
  do 4101 j=1,ndim
 4101 sum=sum+w(i,j)
  do 4102 j=1,ndim-1
 4102 w(i,j)=w(i,j)/sum
  do 4211 j=1,ndim-1
  sum=0.0
  sumsq=0.0
  do 4209 i=1,ns
 4209 sum=sum+w(i,j)
  sum=sum/ns
  do 4210 i=1,ns
  w(i,j)=w(i,j)-sum
 4210 sumsq=sumsq+w(i,j)*w(i,j)
  sumsq=sqrt(sumsq/ns)
  do 4211 i=1,ns
 4211 w(i,j)=w(i,j)/sumsq
  call page ( lout )
  write(lout,4103)title,(j,j=1,ndim-1)
 4103 format(a80/6x,'flattened subject weights'// &
  27x,'variable'/6x,'subject',3x,'plot',i6,5i12)
  write(lout,2683)
  do 4104 i=1,ns
  kk=mod(i-1,52)+1
 4104 write(lout,521)i,alfa(kk),(w(i,j),j=1,ndim-1)
 6566 if ( npt == 0)  go to 9999
  if ( ndim-1 == 1)go to 4108
  do 4107 i=2,ndim-1
  do 4107 j=1,i-1
  call page ( lout )
  write(lout,4106)title,j,i
 4106 format(a80/6x,'flattened subject weights'/6x,'variable', &
  i3,'  (horizontal)  vs  variable',i3,'  (vertical)')
 4107 call plotr(w(1,j),w(1,i),xhi,yhi,xlo,ylo,ns,lout,2,ns)
  go to 9999
 4108 call page ( lout )
  write(lout,4109)title
 4109 format(a80/6x,'flatenned subject weights'/6x, &
  'one variable plot')
  call plotr(w(1,1),w(1,1),xhi,yhi,xlo,ylo,ns,lout,2,ns)
 9999 return
end
subroutine trs ( disp, dist, ivec, nele )
!
!*******************************************************************************
!
!! TRS computes Kruskal's least squares monotonic transformation.
!
!
!  Discussion:
!
!    This routine destroys the distances in DIST, thus forcing
!    double storage in the calling routine.
!
!  Author:
!
!    Rostyslawarema Lewyckyj
!
  dimension disp(1)
  dimension dist(1)
  integer i
  integer ivec(1)
  integer j
  integer nele
!
!  Place distances in ascending order in DIST.
!
  dist(1:nele) = disp(ivec(1:nele))
!
!     perform kruskal's least squares monotonic transformation
!     placing resulting disparities in disp
!
  do i = 2, nele
!
!  Determine if order is correct.
!
  if ( dist(i) >= dist(i-1) ) then
    cycle
  end if

    sum = dist(i)
    fn = 1.0
!
!     if not determine block size
!
    do j = 1, i-1

      imj = i-j
      sum = sum + dist(i-j)
      fn = fn + 1.0
      dispt = sum / fn

      if ( j == i-1 ) then
        exit
      end if

      if ( dispt >= dist(i-j-1) ) then
        exit
      end if

    end do
!
!     set block of disparities equal to mean distance in block
!
    dist(imj:i) = dispt

  end do
!
!  Place disparities in DISP in standard order.
!
  disp(1:nele) = dist(1:nele)

  dist(ivec(1:nele)) = disp(1:nele)

  return
end
subroutine wplot ( cfl, w, ws, nwe, ndim, npoint, ns )
!
!*******************************************************************************
!
!! WPLOT writes file "cplot.dat" for rotation of configuration by use of CPLOT.
!
!
!   wplot is called by step4
!
  character ( len = 80 ) fmt
  character ( len = 80 ) title
!
  common /block3/ title,fmt
  common /cpplot/nrow
!
  dimension cfl(npoint,1),w(ns,1),ws(npoint,1)
  nplt=99

  open ( unit = nplt, file = 'cplot.dat' )

  write (nplt,10) title
   10 format (1x,a)
!
!  default values for parameters
!
  if ( nrow < npoint ) then
    ncol = npoint - nrow
  else
    ncol = 0
  end if

  write (nplt,20) ndim,npoint,ncol
   20 format (3i4,'  1  0 0 0 0',/' data')
!
!  write configuration
!
  do i = 1, npoint
    write(nplt,251) cfl(i,1:ndim)
  end do
  251 format(6f9.5)
!
!  write sample target
!
  i=1
  write (nplt,30) i,(j,j=1,npoint)
   30 format (' target',/,i3,2x,2(25i3))
  write(nplt,251) (cfl(1,j),j=1,ndim)

  if ( nwe/2*2 /= nwe) then
    write(nplt,*)
    do i=1,ns
      write(nplt,251) (w(i,j),j=1,ndim)
    end do
  end if

  if ( nwe >= 2) then
    write(nplt,*)
    do i=1,npoint
      write(nplt,251) (ws(i,j),j=1,ndim)
    end do
  end if

  close ( unit = nplt )
  return
end
