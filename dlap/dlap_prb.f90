program main

!*****************************************************************************80
!
!! MAIN is the main program for DLAP_PRB.
!
!  Discussion:
!
!    DLAP_PRB runs the quick checks for the DLAP sparse linear algebra package.
!
!    This is a SLATEC Quick Checks program to test the *SLAP*
!    Version 2.0 package.  It generates a "random" matrix (See
!    DRMGEN) and then runs all the various methods with all the
!    various preconditoners and all the various stop tests.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager
!    Lawrence Livermore National Laboratory
!
!  Local Parameters:
!
!    Local, integer KPRINT, determines the amount of output.
!    0  Quick checks - No printing.
!       Driver       - Short pass or fail message printed.
!    1  Quick checks - No message printed for passed tests,
!                      short message printed for failed tests.
!       Driver       - Short pass or fail message printed.
!    2  Quick checks - Print short message for passed tests,
!                      fuller information for failed tests.
!       Driver       - Pass or fail message printed.
!    3  Quick checks - Print complete quick check results.
!       Driver       - Pass or fail message printed.
!    4  Quick checks - Print complete quick check results.
!                      Prints matricies, etc.  Very verbose.
!       Driver       - Pass or fail message printed.
!
  implicit double precision (a-h,o-z)

  integer, parameter :: maxiw = 50000
  integer, parameter :: maxn = 441
  integer, parameter :: maxrw = 50000
  integer, parameter :: mxnelt = 50000

  double precision a(mxnelt)
  double precision f(maxn)
  integer ia(mxnelt)
  integer iwork(maxiw)
  integer ja(mxnelt)
  integer, parameter :: kprint = 2
  character ( len = 72 ) mesg
  double precision rwork(maxrw)
  double precision xiter(maxn)
!
!  Supplying the exact solution in common block SOLBLK allows DLAP,
!  when ITOL = 11, to compare its computed answer to the exact solution.
!
  common /solblk/ soln(maxn)

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DLAP_PRB:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the DLAP library.'
!
!  Set up the error routines.
!
  istdi = i1mach(1)
  istdo = i1mach(2)
  nfail = 0

  call xsetun ( lun )

  if ( kprint <= 1 ) then
    call xsetf ( 0 )
  else
    call xsetf ( 1 )
  end if

  call xermax ( 1000 )
!
!  Set the maximum problem sizes.
!
  neltmx = mxnelt
  nmax   = maxn
  leniw  = maxiw
  lenw   = maxrw
!
!  Set some input data.
!
  n      = nmax
  itmax  = n
  iout   = kprint
  factor = 1.2D+00

  if ( iout < 3 ) then
    iunit = 0
  else
    iunit = istdo
  end if
!
!  Set the error tolerance to depend on the machine epsilon.
!
  tol = max ( 1.0D+03 * epsilon ( tol ), 1.0D-06 )
!
!  Test routines using various convergence criteria.
!
  do kase = 3, 3

    if ( kase  ==  1 .or. kase  ==  2 ) then
      itol = kase
    else if ( kase  ==  3 ) then
      itol = 11
    end if
!
!  Test routines using nonsymmetric (ISYM=0) and symmetric
!  storage (ISYM=1).  For ISYM=0 a really non-symmetric matrix
!  is generated.  The amount of non-symmetry is controlled by
!  user.
!
     do isym = 0, 1
!
!  Set up a random matrix.
!
        call drmgen ( neltmx, factor, ierr, n, nelt, &
          isym, ia, ja, a, f, soln, rwork, iwork, iwork(n+1) )

        if ( ierr /= 0 ) then
           mesg = 'DLAPQC -- Fatal error (i1) generating '// &
                '*random* matrix.'
           call xerrwv( mesg,len(mesg),ierr,2,1,ierr,0, &
                0,0.0,0.0 )
        end if

        if ( isym == 0 ) then
           dens = real ( nelt, kind = 8 ) / real ( n*n, kind = 8 )
        else
           dens = real ( 2*nelt, kind = 8 ) / real ( n*n, kind = 8 )
        end if

        if ( 2 <= iout ) then
          write(istdo,1020) n, nelt, dens
          write(istdo,1030) tol
        end if
!
!  Convert to the SLAP-Column format and
!  write out matrix in SLAP-Column format, if desired.
!
        call ds2y( n, nelt, ia, ja, a, isym )

        if ( 4 <= iout ) then
           write(istdo,1040) (k,ia(k),ja(k),a(k),k=1,nelt)
           call dcpplt( n, nelt, ia, ja, a, isym, istdo )
        end if
!
!  BEGIN SLAP QUICK TESTS
!
!  DSJAC.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSJAC :  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsjac(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, 2*itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dsjac ',ierr,iout,nfail,istdo,iter,err )
!
!  DSGS.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSOS  :  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsgs(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dsgs  ',ierr,iout,nfail,istdo,iter,err )
!
!  DSILUR.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSILUR:  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsilur(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dsilur',ierr,iout,nfail,istdo,iter,err )
!
!  DSDCG.
!
        if ( isym == 1 ) then

           if ( 3 <= iout ) then
              write ( *, '(a)' ) ' '
              write ( *, '(2x,a6,a, i2,a,i1 )' ) &
                'DSDCG :  ITOL = ', itol, '  ISYM = ', isym
           end if

           xiter(1:n) = 0.0D+00

           call dsdcg(n, f, xiter, nelt, ia, ja, a, isym, &
                itol, tol, itmax, iter, err, ierr, iunit, &
                rwork, lenw, iwork, leniw )

           call duterr( 'dsdcg ',ierr,iout,nfail,istdo,iter,err )
        end if
!
!  DSICCG.
!
        if ( isym == 1 ) then

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1 )' ) &
               'DSICG :  ITOL = ', itol, '  ISYM = ', isym
           end if

           xiter(1:n) = 0.0D+00

           call dsiccg(n, f, xiter, nelt, ia, ja, a, isym, &
                itol, tol, itmax, iter, err, ierr, iunit, rwork, &
                lenw, iwork, leniw )

           call duterr( 'dsiccg',ierr,iout,nfail,istdo,iter,err )

        end if
!
!  DSDCGN.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSDCGN:  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dsdcgn(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dsdcgn',ierr,iout,nfail,istdo,iter,err )
!
!  DSLUCN.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'DSLUCN:  ITOL = ', itol, '  ISYM = ', isym
        end if

        xiter(1:n) = 0.0D+00

        call dslucn(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dslucn',ierr,iout,nfail,istdo,iter,err )
!
!  DSDBCG.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dsdbcg:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dsdbcg(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dsdbcg',ierr,iout,nfail,istdo,iter,err )
!
!  DSLUBC.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dslubc:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dslubc(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dslubc',ierr,iout,nfail,istdo,iter,err )
!
!  DSDCGS.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dsdcos:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dsdcgs(n, f, xiter, nelt, ia, ja, a, isym, itol, &
             tol, itmax, iter, err, ierr, iunit, rwork, lenw, &
             iwork, leniw )

        call duterr( 'dsdcgs',ierr,iout,nfail,istdo,iter,err )
!
!  DSLUCS.
!
        if ( 3 <= iout ) then
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a6,a, i2,a,i1 )' ) &
            'dslucs:  itol = ', itol, '  isym = ', isym
        end if

        xiter(1:n) = 0.0d+00

        call dslucs(n, f, xiter, nelt, ia, ja, a, isym, &
             itol, tol, itmax, iter, err, ierr, iunit, &
             rwork, lenw, iwork, leniw )

        call duterr( 'dslucs',ierr,iout,nfail,istdo,iter,err )
!
!  DSDOMN.
!
!VD$ NOVECTOR

        do nsave = 0, 3

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dsdomn:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00

           call dsdomn(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itol, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dsdomn',ierr,iout,nfail,istdo,iter,err )

        end do
!
!  DSLUOM
!
!VD$ NOVECTOR

        do nsave = 0, 3

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dsluom:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00

           call dsluom(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itol, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dsluom',ierr,iout,nfail,istdo,iter,err )

        end do
!
!  DSDGMR
!
!VD$ NOVECTOR

        do nsave = 5, 12

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dsdgmr:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00
           itolgm = 0

           call dsdgmr(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itolgm, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dsdgmr',ierr,iout,nfail,istdo,iter,err )

        end do
!
!  DSLUGM
!
!VD$ NOVECTOR

        do nsave = 5, 12

           if ( 3 <= iout ) then
             write ( *, '(a)' ) ' '
             write ( *, '(2x,a6,a, i2,a,i1,a,i2 )' ) &
               'dslugm:  itol = ', itol, '  isym = ', isym, &
               '  nsave = ', nsave
           end if

           xiter(1:n) = 0.0d+00

           call dslugm(n, f, xiter, nelt, ia, ja, a, &
                isym, nsave, itol, tol, itmax, iter, err, ierr, &
                iunit, rwork, lenw, iwork, leniw )

           call duterr( 'dslugm',ierr,iout,nfail,istdo,iter,err )

        end do

     end do

  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '

  if ( NFAIL == 0 ) then
    write ( *, '(a)' ) '*******************************************************'
    write ( *, '(a)' ) '**** All SLAP Double Precision Quick Checks Passed ****'
    write ( *, '(a)' ) '****                 No Errors                     ****'
    write ( *, '(a)' ) '*******************************************************'
  else
     write(istdo,1060) nfail
  end if
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DLAP_PRB:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop

 1020 FORMAT(/'                * RANDOM Matrix of size',I5,'*' &
       /'                ', &
       'Number of non-zeros & Density = ', I5,E16.7)
 1030 FORMAT('                Error tolerance = ',E16.7)
 1040 FORMAT(/'  ***** SLAP Column Matrix *****'/ &
          ' Indx   ia   ja     a'/(1X,I4,1X,I4,1X,I4,1X,E16.7))

 1060 FORMAT(// &
       '************************************************'/ &
       '**     ===>',I3,' Failures detected <===      **'/ &
       '**     SLAP Double Precision Quick Checks     **'/ &
       '** Set KPRINT = 3 for DEBUG information and   **'/ &
       '** rerun the tests to determine the problem.  **'/ &
       '************************************************')
end
subroutine DUTERR ( METHOD, IERR, IOUT, NFAIL, ISTDO, ITER, ERR )

!*****************************************************************************80
!
!! DUTERR outputs error messages for the SLAP quick checks.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager,
!    Lawrence Livermore National Laboratory
!
!  Parameters:
!
  implicit double precision(a-h,o-z)

  double precision err
  integer ierr
  integer iout
  integer istdo
  integer iter
  character*6 method
  integer nfail

  if ( ierr /= 0 ) then
    nfail = nfail+1
  end if

  if ( iout == 1 .and. ierr /= 0 ) then
     write(istdo,1000) method
  end if

  if ( iout == 2 ) then
     if ( ierr == 0 ) then
        write(istdo,1010) method
     else
        write(istdo,1020) method,ierr,iter,err
     end if
  end if

  if ( 3 <= iout ) then
     if ( ierr == 0 ) then
        write(istdo,1030) method,ierr,iter,err
     else
        write(istdo,1020) method,ierr,iter,err
     end if
  end if

  return
 1000 FORMAT( 1X,A6,' : **** FAILURE ****')
 1010 FORMAT( 1X,A6,' : **** PASSED  ****')
 1020 FORMAT(' **************** WARNING ***********************'/ &
         ' **** ',A6,' Quick Test FAILED: IERR = ',I5,' ****'/ &
         ' **************** WARNING ***********************'/ &
         ' Iteration Count = ',I3,' Stop Test = ',E12.6)
 1030 FORMAT(' ***************** PASSED ***********************'/ &
         ' **** ',A6,' Quick Test PASSED: IERR = ',I5,' ****'/ &
         ' ***************** PASSED ***********************'/ &
         ' Iteration Count = ',I3,' Stop Test = ',E12.6)
end
subroutine drmgen ( neltmx, factor, ierr, n, nelt, isym, &
  ia, ja, a, f, soln, dsum, itmp, idiag )

!*****************************************************************************80
!
!! DRMGEN generates a random matrix for SLAP quick checks.
!
!  Discussion:
!
!    The matrix is generated by choosing a random number of
!    entries for each column and then chosing negative random
!    numbers for each off diagionals.   The diagionals elements
!    are chosen to be positive and large enough so the matrix
!    is slightly diagionally dominant.  The lower triangle of
!    the matrix is generated and if isym == 0 (all matrix elements
!    stored) the upper triangle elements are chosen so that they
!    are FACTOR times the coresponding lower triangular element.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager,
!    Lawrence Livermore National Laboratory
!
!  Parameters:
!
!    NELTMX :IN       Integer.
!    Maximum number of non-zeros that can be created by this
!    routine for storage in the IA, JA, A arrays,  see below.
!
!    FACTOR :IN       Double Precision.
!    Non-zeros in the upper triangle are set to FACTOR times
!    the coresponding entry in the lower triangle when a non-
!    symmetric matrix is requested (See ISYM, below).
!
!    IERR   :OUT      Integer.
!    Return error flag.
!    0 => everything went OK.
!    1 => Ran out of space trying to create matrix.
!    Set NELTMX to something larger and retry.
!
!    N      :IN       Integer.
!    Size of the linear system to generate (number of unknowns).
!
!    NELT   :OUT      Integer.
!    Number of non-zeros stored in the IA, JA, A arrays, see below.
!
!    ISYM   :IN       Integer.
!    Flag to indicate the type of matrix to generate:
!    0 => Non-Symmetric Matrix (See FACTOR, above).
!    1 => Symmetric Matrix.
!
!    IA     :OUT      Integer IA(NELTMX).
!    Stores the row indicies for the non-zeros.
!
!    JA     :OUT      Integer JA(NELTMX).
!    Stores the column indicies for the non-zeros.
!
!    A      :OUT      Double Precision A(NELTMX).
!    Stores the values of the non-zeros.
!
!    F      :OUT      Double Precision F(N).
!    The right hand side of the linear system.  Obtained by mult-
!    iplying the matrix time SOLN, see below.
!
!    SOLN   :OUT      Double Precision SOLN(N).
!    The true solution to the linear system.  Each component is
!    chosen at random (0.0<SOLN(I)<1.0, I=1,N)
!
!    DSUM   :WORK     Double Precision DSUM(N).
!
!    ITMP   :WORK     Integer ITMP(N).
!
!    IDIAG  :WORK     Integer IDIAG(N).
!
  implicit double precision(a-h,o-z)

  integer n
  integer neltmx

  double precision a(neltmx)
  double precision dsum(n)
  double precision f(n)
  double precision factor
  integer ia(neltmx)
  integer icol
  integer idiag(n)
  integer ierr
  integer inum
  integer iseed
  integer isym
  integer itmp(n)
  integer ja(neltmx)
  integer nelt
  integer nl
  real rand
  real rn
  double precision soln(n)
!
!  Start by setting the random number generator seed.
!  This is done for reproducability in debugging.  Remove
!  the seed setting call for production testing.
!
  ierr = 0

  idiag(1:n) = 0
  dsum(1:n) = -1.0d+00
!
!  Set the matrix elements.
!  Loop over the columns.
!
  nelt = 0

!VD$ NOCONCUR

  do icol = 1, n

     nl = n+1-icol
!
!  To keep things sparse divide by two, three or four or ...
!
     call random_number ( harvest = rn )

     inum = ( int ( rn * nl ) + 1)/3

     call dmpl ( nl, inum, itmp )
!
!  Set up this column (and row, if non-symmetric structure).
!
!VD$ NOVECTOR
!VD$ NOCONCUR

     do irow = 1, inum

        nelt = nelt + 1

        if ( nelt > neltmx ) then
           ierr = 1
           return
        end if

        ia(nelt) = n+1-itmp(irow)
        ja(nelt) = icol

        if ( ia(nelt) == icol ) then

           idiag(icol) = nelt
        else

           call random_number ( harvest = rn )
           a(nelt) = - rn
           dsum(icol) = dsum(icol) + a(nelt)
           if ( isym == 0 ) then
!
!  Copy this element into upper triangle.
!
              nelt = nelt + 1

              if ( nelt > neltmx ) then
                 ierr = 1
                 return
              end if

              ia(nelt) = icol
              ja(nelt) = ia(nelt-1)
              a(nelt)  = a(nelt-1)*factor
              dsum(ja(nelt)) = dsum(ja(nelt)) + a(nelt)
           else
              dsum(ia(nelt)) = dsum(ia(nelt)) + a(nelt)
           end if

        end if
     end do

     if ( IDIAG(ICOL) == 0 ) then
!
!  Add a diagional to the column.
!
        nelt = nelt + 1
        if ( nelt > neltmx ) then
           ierr = 1
           return
        end if

        idiag(icol) = nelt
        a(nelt) = 0.0d0
        ia(nelt) = icol
        ja(nelt) = icol

     end if

  end do
!
!  Clean up the diagionals.
!
!VD$ NODEPCHK
!LLL. OPTION ASSERT (NOHAZARD)
!DIR$ IVDEP

  do i = 1, n
     a(idiag(i)) = -1.0001D+00 * dsum(i)
  end do
!
!  Set a random solution and determine the right-hand side.
!
!VD$ NOVECTOR
!VD$ NOCONCUR

  call random_number ( harvest = soln(1:n) )

  f(1:n) = 0.0d+00
!
!VD$ NOVECTOR
!VD$ NOCONCUR

  do k = 1, nelt
     f(ia(k)) = f(ia(k)) + a(k)*soln(ja(k))
     if ( isym /= 0 .and. ia(k) /= ja(k) ) then
        f(ja(k)) = f(ja(k)) + a(k)*soln(ia(k))
     end if
  end do

  return
end
subroutine dmpl ( n, m, indx )

!*****************************************************************************80
!
!! DMPL picks M distinct integers at random between 1 and N.
!
!  Modified:
!
!    08 August 2006
!
!  Author:
!
!    Mark Seager
!
!  Parameters:
!
!    Input, integer N, the upper limit on the values.
!
!    Input, integer M, the number of values to pick.
!
!    Output, integer INDX(M), distinct integers between 1 and N.
!
  implicit none

  integer m

  logical found
  integer i
  integer id
  integer indx(m)
  integer j
  integer n
  real rn
!
!  Check the input.
!
  if ( n * m < 0 .or. n < m ) then
    return
  end if
!
!  Set the indices.
!
  call random_number ( harvest = rn )

  indx(1) = int ( rn * real ( n, kind = 8 ) ) + 1

!VD$ NOCONCUR

  do i = 2, m

    do

      call random_number ( harvest = rn )
      id = int ( rn * real ( n, kind = 8 ) ) + 1
!
!  Check to see if id has already been chosen.
!
!VD$ NOVECTOR
!VD$ NOCONCUR

      found = .true.

      do j = 1, i-1

        if ( id == indx(j) ) then
          found = .false.
          exit
        end if

      end do

      if ( found ) then
        indx(i) = id
        exit
      end if

    end do

  end do

  return
end
