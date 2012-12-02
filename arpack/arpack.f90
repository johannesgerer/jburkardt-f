      subroutine sgetv0 &
         ( ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, &
           ipntr, workd, ierr )
!
!! SGETV0 generates a random initial residual vector for the Arnoldi process.
!
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: sgetv0
!
!\Description:
!  SGETV0 generates a random initial residual vector for the Arnoldi process.
!  Force the residual vector to be in the range of the operator OP.
!
!\Usage:
!  call sgetv0
!     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
!       IPNTR, WORKD, IERR )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to sgetv0.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B in the (generalized)
!          eigenvalue problem A*x = lambda*B*x.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  ITRY    Integer.  (INPUT)
!          ITRY counts the number of times that sgetv0 is called.
!          It should be set to 1 on the initial call to sgetv0.
!
!  INITV   Logical variable.  (INPUT)
!          .TRUE.  => the initial residual vector is given in RESID.
!          .FALSE. => generate a random initial residual vector.
!
!  N       Integer.  (INPUT)
!          Dimension of the problem.
!
!  J       Integer.  (INPUT)
!          Index of the residual vector to be generated, with respect to
!          the Arnoldi process.  J > 1 in case of a "restart".
!
!  V       Real N by J array.  (INPUT)
!          The first J-1 columns of V contain the current Arnoldi basis
!          if this is a "restart".
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          Initial residual vector to be generated.  If RESID is
!          provided, force RESID into the range of the operator OP.
!
!  RNORM   Real scalar.  (OUTPUT)
!          B-norm of the generated residual.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!
!  WORKD   Real work array of length 2*N.  (REVERSE COMMUNICATION).
!          On exit, WORK(1:N) = B*RESID to be used in SSAITR.
!
!  IERR    Integer.  (OUTPUT)
!          =  0: Normal exit.
!          = -1: Cannot generate a nontrivial restarted residual vector
!                in the range of the operator OP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     svout   ARPACK utility routine for vector output.
!     slarnv  LAPACK routine for generating a random vector.
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     scopy   Level 1 BLAS that copies one vector to another.
!     sdot    Level 1 BLAS that computes the scalar product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: getv0.F   SID: 2.7   DATE OF SID: 04/07/99   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1
      logical    initv
      integer    ido, ierr, itry, j, ldv, n
      Real &
                 rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Real &
                 resid(n), v(ldv,j), workd(2*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    first, inits, orth
      integer    idist, iseed(4), iter, msglvl, jj
      Real &
                 rnorm0
      save       first, iseed, inits, iter, msglvl, orth, rnorm0
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   slarnv, svout, scopy, sgemv
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 sdot, snrm2
      external   sdot, snrm2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %-----------------%
!     | Data Statements |
!     %-----------------%
!
      data       inits /.true./
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!
!     %-----------------------------------%
!     | Initialize the seed of the LAPACK |
!     | random number generator           |
!     %-----------------------------------%
!
      if (inits) then
          iseed(1) = 1
          iseed(2) = 3
          iseed(3) = 5
          iseed(4) = 7
          inits = .false.
      end if
!
      if (ido ==  0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call cpu_time (t0)
         msglvl = mgetv0
!
         ierr   = 0
         iter   = 0
         first  = .FALSE.
         orth   = .FALSE.
!
!        %-----------------------------------------------------%
!        | Possibly generate a random starting vector in RESID |
!        | Use a LAPACK random number generator used by the    |
!        | matrix generation routines.                         |
!        |    idist = 1: uniform (0,1)  distribution;          |
!        |    idist = 2: uniform (-1,1) distribution;          |
!        |    idist = 3: normal  (0,1)  distribution;          |
!        %-----------------------------------------------------%
!
         if (.not.initv) then
            idist = 2
            call slarnv (idist, iseed, n, resid)
         end if
!
!        %----------------------------------------------------------%
!        | Force the starting vector into the range of OP to handle |
!        | the generalized problem when B is possibly (singular).   |
!        %----------------------------------------------------------%
!
         call cpu_time (t2)
         if (bmat == 'G') then
            nopx = nopx + 1
            ipntr(1) = 1
            ipntr(2) = n + 1
            call scopy (n, resid, 1, workd, 1)
            ido = -1
            go to 9000
         end if
      end if
!
!     %-----------------------------------------%
!     | Back from computing OP*(initial-vector) |
!     %-----------------------------------------%
!
      if (first) go to 20
!
!     %-----------------------------------------------%
!     | Back from computing B*(orthogonalized-vector) |
!     %-----------------------------------------------%
!
      if (orth)  go to 40
!
      if (bmat == 'G') then
         call cpu_time (t3)
         tmvopx = tmvopx + (t3 - t2)
      end if
!
!     %------------------------------------------------------%
!     | Starting vector is now in the range of OP; r = OP*r; |
!     | Compute B-norm of starting vector.                   |
!     %------------------------------------------------------%
!
      call cpu_time (t2)
      first = .TRUE.
      if (bmat == 'G') then
         nbx = nbx + 1
         call scopy (n, workd(n+1), 1, resid, 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat == 'I') then
         call scopy (n, resid, 1, workd, 1)
      end if
!
   20 continue
!
      if (bmat == 'G') then
         call cpu_time (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
!
      first = .FALSE.
      if (bmat == 'G') then
          rnorm0 = sdot (n, resid, 1, workd, 1)
          rnorm0 = sqrt(abs(rnorm0))
      else if (bmat == 'I') then
           rnorm0 = snrm2(n, resid, 1)
      end if
      rnorm  = rnorm0
!
!     %---------------------------------------------%
!     | Exit if this is the very first Arnoldi step |
!     %---------------------------------------------%
!
      if (j == 1) go to 50
!
!     %----------------------------------------------------------------
!     | Otherwise need to B-orthogonalize the starting vector against |
!     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
!     | This is the case where an invariant subspace is encountered   |
!     | in the middle of the Arnoldi factorization.                   |
!     |                                                               |
!     |       s = V^{T}*B*r;   r = r - V*s;                           |
!     |                                                               |
!     | Stopping criteria used for iter. ref. is discussed in         |
!     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
!     %---------------------------------------------------------------%
!
      orth = .TRUE.
   30 continue
!
      call sgemv ('T', n, j-1, one, v, ldv, workd, 1, &
                  zero, workd(n+1), 1)
      call sgemv ('N', n, j-1, -one, v, ldv, workd(n+1), 1, &
                  one, resid, 1)
!
!     %----------------------------------------------------------%
!     | Compute the B-norm of the orthogonalized starting vector |
!     %----------------------------------------------------------%
!
      call cpu_time (t2)
      if (bmat == 'G') then
         nbx = nbx + 1
         call scopy (n, resid, 1, workd(n+1), 1)
         ipntr(1) = n + 1
         ipntr(2) = 1
         ido = 2
         go to 9000
      else if (bmat == 'I') then
         call scopy (n, resid, 1, workd, 1)
      end if
!
   40 continue
!
      if (bmat == 'G') then
         call cpu_time (t3)
         tmvbx = tmvbx + (t3 - t2)
      end if
!
      if (bmat == 'G') then
         rnorm = sdot (n, resid, 1, workd, 1)
         rnorm = sqrt(abs(rnorm))
      else if (bmat == 'I') then
         rnorm = snrm2(n, resid, 1)
      end if
!
!     %--------------------------------------%
!     | Check for further orthogonalization. |
!     %--------------------------------------%
!
      if (msglvl > 2) then
          call svout (logfil, 1, rnorm0, ndigit, &
                      '_getv0: re-orthonalization ; rnorm0 is')
          call svout (logfil, 1, rnorm, ndigit, &
                      '_getv0: re-orthonalization ; rnorm is')
      end if
!
      if (rnorm > 0.717*rnorm0) go to 50
!
      iter = iter + 1
      if (iter <= 5) then
!
!        %-----------------------------------%
!        | Perform iterative refinement step |
!        %-----------------------------------%
!
         rnorm0 = rnorm
         go to 30
      else
!
!        %------------------------------------%
!        | Iterative refinement step "failed" |
!        %------------------------------------%
!
         do 45 jj = 1, n
            resid(jj) = zero
   45    continue
         rnorm = zero
         ierr = -1
      end if
!
   50 continue
!
      if (msglvl > 0) then
         call svout (logfil, 1, rnorm, ndigit, &
              '_getv0: B-norm of initial / restarted starting vector')
      end if
      if (msglvl > 3) then
         call svout (logfil, n, resid, ndigit, &
              '_getv0: initial / restarted starting vector')
      end if
      ido = 99
!
      call cpu_time (t1)
      tgetv0 = tgetv0 + (t1 - t0)
!
 9000 continue
      return
!
!     %---------------%
!     | End of sgetv0 |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: slaqrb
!
!\Description:
!  SLAQRB computes eigenvalues and Schur decomposition of an upper
!  Hessenberg submatrix in rows and columns ILO to IHI.  Only the
!  last component of the Schur vectors are computed.
!
!  This is mostly a modification of the LAPACK routine slahqr.
!
!\Usage:
!  call slaqrb
!     ( WANTT, N, ILO, IHI, H, LDH, WR, WI,  Z, INFO )
!
!\Arguments
!  WANTT   Logical variable.  (INPUT)
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!  N       Integer.  (INPUT)
!          The order of the matrix H.  N >= 0.
!
!  ILO     Integer.  (INPUT)
!  IHI     Integer.  (INPUT)
!          It is assumed that H is already upper quasi-triangular in
!          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!          ILO = 1). SLAQRB works primarily with the Hessenberg
!          submatrix in rows and columns ILO to IHI, but applies
!          transformations to all of H if WANTT is .TRUE..
!          1 <= ILO <= max(1,IHI); IHI <= N.
!
!  H       Real array, dimension (LDH,N).  (INPUT/OUTPUT)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
!          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
!          standard form. If WANTT is .FALSE., the contents of H are
!          unspecified on exit.
!
!  LDH     Integer.  (INPUT)
!          The leading dimension of the array H. LDH >= max(1,N).
!
!  WR      Real array, dimension (N).  (OUTPUT)
!  WI      Real array, dimension (N).  (OUTPUT)
!          The real and imaginary parts, respectively, of the computed
!          eigenvalues ILO to IHI are stored in the corresponding
!          elements of WR and WI. If two eigenvalues are computed as a
!          complex conjugate pair, they are stored in consecutive
!          elements of WR and WI, say the i-th and (i+1)th, with
!          WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
!          eigenvalues are stored in the same order as on the diagonal
!          of the Schur form returned in H, with WR(i) = H(i,i), and, if
!          H(i:i+1,i:i+1) is a 2-by-2 diagonal block,
!          WI(i) = sqrt(H(i+1,i)*H(i,i+1)) and WI(i+1) = -WI(i).
!
!  Z       Real array, dimension (N).  (OUTPUT)
!          On exit Z contains the last components of the Schur vectors.
!
!  INFO    Integer.  (OUPUT)
!          = 0: successful exit
!          > 0: SLAQRB failed to compute all the eigenvalues ILO to IHI
!               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
!               elements i+1:ihi of WR and WI contain those eigenvalues
!               which have been successfully computed.
!
!\Remarks
!  1. None.
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     slabad  LAPACK routine that computes machine constants.
!     slamch  LAPACK routine that determines machine constants.
!     slanhs  LAPACK routine that computes various norms of a matrix.
!     slanv2  LAPACK routine that computes the Schur factorization of
!             2 by 2 nonsymmetric matrix in standard form.
!     slarfg  LAPACK Householder reflection construction routine.
!     scopy   Level 1 BLAS that copies one vector to another.
!     srot    Level 1 BLAS that applies a rotation to a 2 by 2 matrix.

!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!               Modified from the LAPACK routine slahqr so that only the
!               last component of the Schur vectors are computed.
!
!\SCCS Information: @(#)
! FILE: laqrb.F   SID: 2.2   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine slaqrb ( wantt, n, ilo, ihi, h, ldh, wr, wi, &
                          z, info )
!
!! SLAQRB computes eigenvalues and Schur decomposition.
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      logical    wantt
      integer    ihi, ilo, info, ldh, n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 h( ldh, * ), wi( * ), wr( * ), z( * )
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 zero, one, dat1, dat2
      parameter (zero = 0.0E+0, one = 1.0E+0, dat1 = 7.5E-1, &
                 dat2 = -4.375E-1)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      integer    i, i1, i2, itn, its, j, k, l, m, nh, nr
      Real &
                 cs, h00, h10, h11, h12, h21, h22, h33, h33s, &
                 h43h34, h44, h44s, ovfl, s, smlnum, sn, sum, &
                 t1, t2, t3, tst1, ulp, unfl, v1, v2, v3
      Real &
                 v( 3 ), work( 1 )
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slamch, slanhs
      external   slamch, slanhs
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy, slabad, slanv2, slarfg, srot
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      info = 0
!
!     %--------------------------%
!     | Quick return if possible |
!     %--------------------------%
!
      if( n==0 ) &
         return
      if( ilo==ihi ) then
         wr( ilo ) = h( ilo, ilo )
         wi( ilo ) = zero
         return
      end if
!
!     %---------------------------------------------%
!     | Initialize the vector of last components of |
!     | the Schur vectors for accumulation.         |
!     %---------------------------------------------%
!
      z(1:n-1) = zero
      z(n) = one
!
      nh = ihi - ilo + 1
!
!     %-------------------------------------------------------------%
!     | Set machine-dependent constants for the stopping criterion. |
!     | If norm(H) <= sqrt(OVFL), overflow should not occur.        |
!     %-------------------------------------------------------------%
!
      unfl = slamch( 'safe minimum' )
      ovfl = one / unfl
      call slabad( unfl, ovfl )
      ulp = slamch( 'precision' )
      smlnum = unfl*( nh / ulp )
!
!     %---------------------------------------------------------------%
!     | I1 and I2 are the indices of the first row and last column    |
!     | of H to which transformations must be applied. If eigenvalues |
!     | only are computed, I1 and I2 are set inside the main loop.    |
!     | Zero out H(J+2,J) = ZERO for J=1:N if WANTT = .TRUE.          |
!     | else H(J+2,J) for J=ILO:IHI-ILO-1 if WANTT = .FALSE.          |
!     %---------------------------------------------------------------%
!
      if( wantt ) then
         i1 = 1
         i2 = n
         do 8 i=1,i2-2
            h(i1+i+1,i) = zero
 8       continue
      else
         do 9 i=1, ihi-ilo-1
            h(ilo+i+1,ilo+i-1) = zero
 9       continue
      end if
!
!     %---------------------------------------------------%
!     | ITN is the total number of QR iterations allowed. |
!     %---------------------------------------------------%
!
      itn = 30*nh
!
!     ------------------------------------------------------------------
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!     ------------------------------------------------------------------
!
      i = ihi
   10 continue
      l = ilo
      if( i<ilo ) &
         go to 150

!     %--------------------------------------------------------------%
!     | Perform QR iterations on rows and columns ILO to I until a   |
!     | submatrix of order 1 or 2 splits off at the bottom because a |
!     | subdiagonal element has become negligible.                   |
!     %--------------------------------------------------------------%

      do 130 its = 0, itn
!
!        %----------------------------------------------%
!        | Look for a single small subdiagonal element. |
!        %----------------------------------------------%
!
         do k = i, l + 1, -1
            tst1 = abs( h( k-1, k-1 ) ) + abs( h( k, k ) )
            if( tst1==zero ) &
               tst1 = slanhs( '1', i-l+1, h( l, l ), ldh, work )
            if( abs( h( k, k-1 ) )<=max( ulp*tst1, smlnum ) ) &
               go to 30
         end do

   30    continue
         l = k
         if( l>ilo ) then
!
!           %------------------------%
!           | H(L,L-1) is negligible |
!           %------------------------%
!
            h( l, l-1 ) = zero
         end if
!
!        %-------------------------------------------------------------%
!        | Exit from loop if a submatrix of order 1 or 2 has split off |
!        %-------------------------------------------------------------%
!
         if( l>=i-1 ) &
            go to 140
!
!        %---------------------------------------------------------%
!        | Now the active submatrix is in rows and columns L to I. |
!        | If eigenvalues only are being computed, only the active |
!        | submatrix need be transformed.                          |
!        %---------------------------------------------------------%
!
         if( .not.wantt ) then
            i1 = l
            i2 = i
         end if
!
         if( its==10 .or. its.eq.20 ) then
!
!           %-------------------%
!           | Exceptional shift |
!           %-------------------%
!
            s = abs( h( i, i-1 ) ) + abs( h( i-1, i-2 ) )
            h44 = dat1*s
            h33 = h44
            h43h34 = dat2*s*s
!
         else
!
!           %-----------------------------------------%
!           | Prepare to use Wilkinson's double shift |
!           %-----------------------------------------%
!
            h44 = h( i, i )
            h33 = h( i-1, i-1 )
            h43h34 = h( i, i-1 )*h( i-1, i )
         end if
!
!        %-----------------------------------------------------%
!        | Look for two consecutive small subdiagonal elements |
!        %-----------------------------------------------------%
!
         do 40 m = i - 2, l, -1
!
!           %---------------------------------------------------------%
!           | Determine the effect of starting the double-shift QR    |
!           | iteration at row M, and see if this would make H(M,M-1) |
!           | negligible.                                             |
!           %---------------------------------------------------------%
!
            h11 = h( m, m )
            h22 = h( m+1, m+1 )
            h21 = h( m+1, m )
            h12 = h( m, m+1 )
            h44s = h44 - h11
            h33s = h33 - h11
            v1 = ( h33s*h44s-h43h34 ) / h21 + h12
            v2 = h22 - h11 - h33s - h44s
            v3 = h( m+2, m+1 )
            s = abs( v1 ) + abs( v2 ) + abs( v3 )
            v1 = v1 / s
            v2 = v2 / s
            v3 = v3 / s
            v( 1 ) = v1
            v( 2 ) = v2
            v( 3 ) = v3
            if( m==l ) &
               go to 50
            h00 = h( m-1, m-1 )
            h10 = h( m, m-1 )
            tst1 = abs( v1 )*( abs( h00 )+abs( h11 )+abs( h22 ) )
            if( abs( h10 )*( abs( v2 )+abs( v3 ) )<=ulp*tst1 ) &
               go to 50
   40    continue
   50    continue
!
!        %----------------------%
!        | Double-shift QR step |
!        %----------------------%
!
         do 120 k = m, i - 1
!
!           ------------------------------------------------------------
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!           ------------------------------------------------------------
!
            nr = min( 3, i-k+1 )
            if( k>m ) &
               call scopy( nr, h( k, k-1 ), 1, v, 1 )
            call slarfg( nr, v( 1 ), v( 2 ), 1, t1 )
            if( k>m ) then
               h( k, k-1 ) = v( 1 )
               h( k+1, k-1 ) = zero
               if( k<i-1 ) &
                  h( k+2, k-1 ) = zero
            else if( m>l ) then
               h( k, k-1 ) = -h( k, k-1 )
            end if
            v2 = v( 2 )
            t2 = t1*v2
            if( nr==3 ) then
               v3 = v( 3 )
               t3 = t1*v3
!
!              %------------------------------------------------%
!              | Apply G from the left to transform the rows of |
!              | the matrix in columns K to I2.                 |
!              %------------------------------------------------%
!
               do 60 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j ) + v3*h( k+2, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
                  h( k+2, j ) = h( k+2, j ) - sum*t3
   60          continue
!
!              %----------------------------------------------------%
!              | Apply G from the right to transform the columns of |
!              | the matrix in rows I1 to min(K+3,I).               |
!              %----------------------------------------------------%
!
               do 70 j = i1, min( k+3, i )
                  sum = h( j, k ) + v2*h( j, k+1 ) + v3*h( j, k+2 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
                  h( j, k+2 ) = h( j, k+2 ) - sum*t3
   70          continue
!
!              %----------------------------------%
!              | Accumulate transformations for Z |
!              %----------------------------------%
!
               sum      = z( k ) + v2*z( k+1 ) + v3*z( k+2 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
               z( k+2 ) = z( k+2 ) - sum*t3

            else if( nr==2 ) then
!
!              %------------------------------------------------%
!              | Apply G from the left to transform the rows of |
!              | the matrix in columns K to I2.                 |
!              %------------------------------------------------%
!
               do 90 j = k, i2
                  sum = h( k, j ) + v2*h( k+1, j )
                  h( k, j ) = h( k, j ) - sum*t1
                  h( k+1, j ) = h( k+1, j ) - sum*t2
   90          continue
!
!              %----------------------------------------------------%
!              | Apply G from the right to transform the columns of |
!              | the matrix in rows I1 to min(K+3,I).               |
!              %----------------------------------------------------%
!
               do j = i1, i
                  sum = h( j, k ) + v2*h( j, k+1 )
                  h( j, k ) = h( j, k ) - sum*t1
                  h( j, k+1 ) = h( j, k+1 ) - sum*t2
               end do
!
!              %----------------------------------%
!              | Accumulate transformations for Z |
!              %----------------------------------%
!
               sum      = z( k ) + v2*z( k+1 )
               z( k )   = z( k ) - sum*t1
               z( k+1 ) = z( k+1 ) - sum*t2
            end if
  120    continue

  130 continue
!
!     %-------------------------------------------------------%
!     | Failure to converge in remaining number of iterations |
!     %-------------------------------------------------------%
!
      info = i
      return

  140 continue

      if( l==i ) then
!
!        %------------------------------------------------------%
!        | H(I,I-1) is negligible: one eigenvalue has converged |
!        %------------------------------------------------------%
!
         wr( i ) = h( i, i )
         wi( i ) = zero

      else if( l==i-1 ) then
!
!        %--------------------------------------------------------%
!        | H(I-1,I-2) is negligible;                              |
!        | a pair of eigenvalues have converged.                  |
!        |                                                        |
!        | Transform the 2-by-2 submatrix to standard Schur form, |
!        | and compute and store the eigenvalues.                 |
!        %--------------------------------------------------------%
!
         call slanv2( h( i-1, i-1 ), h( i-1, i ), h( i, i-1 ), &
                      h( i, i ), wr( i-1 ), wi( i-1 ), wr( i ), wi( i ), &
                      cs, sn )

         if( wantt ) then
!
!           %-----------------------------------------------------%
!           | Apply the transformation to the rest of H and to Z, |
!           | as required.                                        |
!           %-----------------------------------------------------%
!
            if( i2>i ) &
               call srot( i2-i, h( i-1, i+1 ), ldh, h( i, i+1 ), ldh, &
                          cs, sn )
            call srot( i-i1-1, h( i1, i-1 ), 1, h( i1, i ), 1, cs, sn )
            sum      = cs*z( i-1 ) + sn*z( i )
            z( i )   = cs*z( i )   - sn*z( i-1 )
            z( i-1 ) = sum
         end if
      end if
!
!     %---------------------------------------------------------%
!     | Decrement number of remaining iterations, and return to |
!     | start of the main loop with new value of I.             |
!     %---------------------------------------------------------%
!
      itn = itn - its
      i = l - 1
      go to 10

  150 continue
      return
!
!     %---------------%
!     | End of slaqrb |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: snaitr
!
!\Description:
!  SNAITR is a reverse communication interface for applying NP additional 
!  steps to a K step nonsymmetric Arnoldi factorization.
!
!  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
!
!          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
!
!  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
!
!          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
!
!  where OP and B are as in snaupd.  The B-norm of r_{k+p} is also
!  computed and returned.
!
!\Usage:
!  call snaitr
!     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH,
!       IPNTR, WORKD, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!                    This is for the restart phase to force the new
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y,
!                    IPNTR(3) is the pointer into WORK for B * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!          When the routine is used in the "shift-and-invert" mode, the
!          vector B * Q is already available and do not need to be
!          recompute in forming OP * Q.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.  See snaupd.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  K       Integer.  (INPUT)
!          Current size of V and H.
!
!  NP      Integer.  (INPUT)
!          Number of additional Arnoldi steps to take.
!
!  NB      Integer.  (INPUT)
!          Blocksize to be used in the recurrence.
!          Only work for NB = 1 right now.  The goal is to have a
!          program that implement both the block and non-block method.
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          On INPUT:  RESID contains the residual vector r_{k}.
!          On OUTPUT: RESID contains the residual vector r_{k+p}.
!
!  RNORM   Real scalar.  (INPUT/OUTPUT)
!          B-norm of the starting residual on input.
!          B-norm of the updated residual r_{k+p} on output.
!
!  V       Real N by K+NP array.  (INPUT/OUTPUT)
!          On INPUT:  V contains the Arnoldi vectors in the first K
!          columns.
!          On OUTPUT: V contains the new NP Arnoldi vectors in the next
!          NP columns.  The first K columns are unchanged.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
!          H is used to store the generated upper Hessenberg matrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORK for
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!
!  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The calling program should not
!          use WORKD as temporary workspace during the iteration !!!!!!
!          On input, WORKD(1:N) = B*RESID and is used to save some
!          computation at the first step.
!
!  INFO    Integer.  (OUTPUT)
!          = 0: Normal exit.
!          > 0: Size of the spanning invariant subspace of OP found.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     sgetv0  ARPACK routine to generate the initial vector.
!     ivout   ARPACK utility routine that prints integers.
!     smout   ARPACK utility routine that prints matrices
!     svout   ARPACK utility routine that prints vectors.
!     slabad  LAPACK routine that computes machine constants.
!     slamch  LAPACK routine that determines machine constants.
!     slascl  LAPACK routine for careful scaling of a matrix.
!     slanhs  LAPACK routine that computes various norms of a matrix.
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     saxpy   Level 1 BLAS that computes a vector triad.
!     sscal   Level 1 BLAS that scales a vector.
!     scopy   Level 1 BLAS that copies one vector to another .
!     sdot    Level 1 BLAS that computes the scalar product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: naitr.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!  The algorithm implemented is:
!
!  restart = .false.
!  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
!  r_{k} contains the initial residual vector even for k = 0;
!  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
!  computed by the calling program.
!
!  betaj = rnorm ; p_{k+1} = B*r_{k} ;
!  For  j = k+1, ..., k+np  Do
!     1) if ( betaj < tol ) stop or restart depending on j.
!        ( At present tol is zero )
!        if ( restart ) generate a new starting vector.
!     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
!        p_{j} = p_{j}/betaj
!     3) r_{j} = OP*v_{j} where OP is defined as in snaupd
!        For shift-invert mode p_{j} = B*v_{j} is already available.
!        wnorm = || OP*v_{j} ||
!     4) Compute the j-th step residual vector.
!        w_{j} =  V_{j}^T * B * OP * v_{j}
!        r_{j} =  OP*v_{j} - V_{j} * w_{j}
!        H(:,j) = w_{j};
!        H(j,j-1) = rnorm
!        rnorm = || r_(j) ||
!        If (rnorm > 0.717*wnorm) accept step and go back to 1)
!     5) Re-orthogonalization step:
!        s = V_{j}'*B*r_{j}
!        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
!        alphaj = alphaj + s_{j};
!     6) Iterative refinement step:
!        If (rnorm1 > 0.717*rnorm) then
!           rnorm = rnorm1
!           accept step and go back to 1)
!        Else
!           rnorm = rnorm1
!           If this is the first time in step 6), go to 5)
!           Else r_{j} lies in the span of V_{j} numerically.
!              Set r_{j} = 0 and rnorm = 0; go to 1)
!        EndIf
!  End Do
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine snaitr &
         (ido, bmat, n, k, np, nb, resid, rnorm, v, ldv, h, ldh, &
          ipntr, workd, info)
!
!! SNAITR applies more steps to a K step nonsymmetric Arnoldi factorization.
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, nb, np
      Real &
                 rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Real &
                 h(ldh,k+np), resid(n), v(ldv,k+np), workd(3*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      logical    first, orth1, orth2, rstart, step3, step4
      integer    ierr, i, infol, ipj, irj, ivj, iter, itry, j, msglvl, &
                 jj
      Real &
                 betaj, ovfl, temp1, rnorm1, smlnum, tst1, ulp, unfl, &
                 wnorm
      save       first, orth1, orth2, rstart, step3, step4, &
                 ierr, ipj, irj, ivj, iter, itry, j, msglvl, ovfl, &
                 betaj, rnorm1, smlnum, ulp, unfl, wnorm
!
!     %-----------------------%
!     | Local Array Arguments |
!     %-----------------------%
!
      Real &
                 xtemp(2)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   saxpy, scopy, sscal, sgemv, sgetv0, slabad, &
                 svout, smout, ivout
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 sdot, snrm2, slanhs, slamch
      external   sdot, snrm2, slanhs, slamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, sqrt
!
!     %-----------------%
!     | Data statements |
!     %-----------------%
!
      data      first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (first) then
!
!        %-----------------------------------------%
!        | Set machine-dependent constants for the |
!        | the splitting and deflation criterion.  |
!        | If norm(H) <= sqrt(OVFL),               |
!        | overflow should not occur.              |
!        | REFERENCE: LAPACK subroutine slahqr     |
!        %-----------------------------------------%
!
         unfl = slamch( 'safe minimum' )
         ovfl = one / unfl
         call slabad( unfl, ovfl )
         ulp = slamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
!
      if (ido == 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call cpu_time (t0)
         msglvl = mnaitr
!
!        %------------------------------%
!        | Initial call to this routine |
!        %------------------------------%
!
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
         j      = k + 1
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
!
!     %-------------------------------------------------%
!     | When in reverse communication mode one of:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | will be .true. when ....                        |
!     | STEP3: return from computing OP*v_{j}.          |
!     | STEP4: return from computing B-norm of OP*v_{j} |
!     | ORTH1: return from computing B-norm of r_{j+1}  |
!     | ORTH2: return from computing B-norm of          |
!     |        correction to the residual vector.       |
!     | RSTART: return from OP computations needed by   |
!     |         sgetv0.                                 |
!     %-------------------------------------------------%
!
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
!
!     %-----------------------------%
!     | Else this is the first step |
!     %-----------------------------%
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%

 1000 continue
!
         if (msglvl > 1) then
            call ivout (logfil, 1, j, ndigit, &
                        '_naitr: generating Arnoldi vector number')
            call svout (logfil, 1, rnorm, ndigit, &
                        '_naitr: B-norm of the current residual is')
         end if
!
!        %---------------------------------------------------%
!        | STEP 1: Check if the B norm of j-th residual      |
!        | vector is zero. Equivalent to determing whether   |
!        | an exact j-step Arnoldi factorization is present. |
!        %---------------------------------------------------%
!
         betaj = rnorm
         if (rnorm > zero) go to 40
!
!           %---------------------------------------------------%
!           | Invariant subspace found, generate a new starting |
!           | vector which is orthogonal to the current Arnoldi |
!           | basis and continue the iteration.                 |
!           %---------------------------------------------------%
!
            if (msglvl > 0) then
               call ivout (logfil, 1, j, ndigit, &
                           '_naitr: ****** RESTART AT STEP ******')
            end if
!
!           %---------------------------------------------%
!           | ITRY is the loop variable that controls the |
!           | maximum amount of times that a restart is   |
!           | attempted. NRSTRT is used by stat.h         |
!           %---------------------------------------------%
!
            betaj  = zero
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
!
!           %--------------------------------------%
!           | If in reverse communication mode and |
!           | RSTART = .true. flow returns here.   |
!           %--------------------------------------%
!
            call sgetv0 (ido, bmat, itry, .false., n, j, v, ldv, &
                         resid, rnorm, ipntr, workd, ierr)
            if (ido /= 99) go to 9000
            if (ierr < 0) then
               itry = itry + 1
               if (itry <= 3) go to 20
!
!              %------------------------------------------------%
!              | Give up after several restart attempts.        |
!              | Set INFO to the size of the invariant subspace |
!              | which spans OP and exit.                       |
!              %------------------------------------------------%
!
               info = j - 1
               call cpu_time (t1)
               tnaitr = tnaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
!
   40    continue
!
!        %---------------------------------------------------------%
!        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
!        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
!        | when reciprocating a small RNORM, test against lower    |
!        | machine bound.                                          |
!        %---------------------------------------------------------%
!
         call scopy (n, resid, 1, v(1,j), 1)
         if (rnorm >= unfl) then
             temp1 = one / rnorm
             call sscal (n, temp1, v(1,j), 1)
             call sscal (n, temp1, workd(ipj), 1)
         else
!
!            %-----------------------------------------%
!            | To scale both v_{j} and p_{j} carefully |
!            | use LAPACK routine SLASCL               |
!            %-----------------------------------------%
!
             call slascl ('General', i, i, rnorm, one, n, 1, &
                          v(1,j), n, infol)
             call slascl ('General', i, i, rnorm, one, n, 1, &
                          workd(ipj), n, infol)
         end if
!
!        %------------------------------------------------------%
!        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
!        | Note that this is not quite yet r_{j}. See STEP 4    |
!        %------------------------------------------------------%
!
         step3 = .true.
         nopx  = nopx + 1
         call cpu_time (t2)
         call scopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
!
!        %-----------------------------------%
!        | Exit in order to compute OP*v_{j} |
!        %-----------------------------------%
!
         go to 9000
   50    continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
!        | if step3 = .true.                |
!        %----------------------------------%
!
         call cpu_time (t3)
         tmvopx = tmvopx + (t3 - t2)

         step3 = .false.
!
!        %------------------------------------------%
!        | Put another copy of OP*v_{j} into RESID. |
!        %------------------------------------------%
!
         call scopy (n, workd(irj), 1, resid, 1)
!
!        %---------------------------------------%
!        | STEP 4:  Finish extending the Arnoldi |
!        |          factorization to length j.   |
!        %---------------------------------------%
!
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %-------------------------------------%
!           | Exit in order to compute B*OP*v_{j} |
!           %-------------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   60    continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
!        | if step4 = .true.                |
!        %----------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         step4 = .false.
!
!        %-------------------------------------%
!        | The following is needed for STEP 5. |
!        | Compute the B-norm of OP*v_{j}.     |
!        %-------------------------------------%
!
         if (bmat == 'G') then
             wnorm = sdot (n, resid, 1, workd(ipj), 1)
             wnorm = sqrt(abs(wnorm))
         else if (bmat == 'I') then
            wnorm = snrm2(n, resid, 1)
         end if
!
!        %-----------------------------------------%
!        | Compute the j-th residual corresponding |
!        | to the j step factorization.            |
!        | Use Classical Gram Schmidt and compute: |
!        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
!        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
!        %-----------------------------------------%
!
!
!        %------------------------------------------%
!        | Compute the j Fourier coefficients w_{j} |
!        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
!        %------------------------------------------%
!
         call sgemv ('T', n, j, one, v, ldv, workd(ipj), 1, &
                     zero, h(1,j), 1)
!
!        %--------------------------------------%
!        | Orthogonalize r_{j} against V_{j}.   |
!        | RESID contains OP*v_{j}. See STEP 3. |
!        %--------------------------------------%
!
         call sgemv ('N', n, j, -one, v, ldv, h(1,j), 1, &
                     one, resid, 1)
!
         if (j > 1) h(j,j-1) = betaj
!
         call cpu_time (t4)
!
         orth1 = .true.
!
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*r_{j} |
!           %----------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH1 = .true. |
!        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
!        %---------------------------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         orth1 = .false.
!
!        %------------------------------%
!        | Compute the B-norm of r_{j}. |
!        %------------------------------%
!
         if (bmat == 'G') then
            rnorm = sdot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat == 'I') then
            rnorm = snrm2(n, resid, 1)
         end if
!
!        %-----------------------------------------------------------%
!        | STEP 5: Re-orthogonalization / Iterative refinement phase |
!        | Maximum NITER_ITREF tries.                                |
!        |                                                           |
!        |          s      = V_{j}^T * B * r_{j}                     |
!        |          r_{j}  = r_{j} - V_{j}*s                         |
!        |          alphaj = alphaj + s_{j}                          |
!        |                                                           |
!        | The stopping criteria used for iterative refinement is    |
!        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
!        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
!        | Determine if we need to correct the residual. The goal is |
!        | to enforce ||v(:,1:j)^T * r_{j}|| <= eps * || r_{j} ||  |
!        | The following test determines whether the sine of the     |
!        | angle between  OP*x and the computed residual is less     |
!        | than or equal to 0.717.                                   |
!        %-----------------------------------------------------------%
!
         if (rnorm > 0.717*wnorm) go to 100
         iter  = 0
         nrorth = nrorth + 1
!
!        %---------------------------------------------------%
!        | Enter the Iterative refinement phase. If further  |
!        | refinement is necessary, loop back here. The loop |
!        | variable is ITER. Perform a step of Classical     |
!        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
!        %---------------------------------------------------%
!
   80    continue
!
         if (msglvl > 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call svout (logfil, 2, xtemp, ndigit, &
                 '_naitr: re-orthonalization; wnorm and rnorm are')
            call svout (logfil, j, h(1,j), ndigit, &
                        '_naitr: j-th column of H')
         end if
!
!        %----------------------------------------------------%
!        | Compute V_{j}^T * B * r_{j}.                       |
!        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
         call sgemv ('T', n, j, one, v, ldv, workd(ipj), 1, &
                     zero, workd(irj), 1)
!
!        %---------------------------------------------%
!        | Compute the correction to the residual:     |
!        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
!        | The correction to H is v(:,1:J)*H(1:J,1:J)  |
!        | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
!        %---------------------------------------------%
!
         call sgemv ('N', n, j, -one, v, ldv, workd(irj), 1, &
                     one, resid, 1)
         call saxpy (j, one, workd(irj), 1, h(1,j), 1)
!
         orth2 = .true.
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %-----------------------------------%
!           | Exit in order to compute B*r_{j}. |
!           | r_{j} is the corrected residual.  |
!           %-----------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH2 = .true. |
!        %---------------------------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
!        %-----------------------------------------------------%
!        | Compute the B-norm of the corrected residual r_{j}. |
!        %-----------------------------------------------------%
!
         if (bmat == 'G') then
             rnorm1 = sdot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat == 'I') then
             rnorm1 = snrm2(n, resid, 1)
         end if

         if (msglvl > 0 .and. iter .gt. 0) then
            call ivout (logfil, 1, j, ndigit, &
                 '_naitr: Iterative refinement for Arnoldi residual')
            if (msglvl > 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call svout (logfil, 2, xtemp, ndigit, &
                 '_naitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
!
!        %-----------------------------------------%
!        | Determine if we need to perform another |
!        | step of re-orthogonalization.           |
!        %-----------------------------------------%
!
         if (rnorm1 > 0.717*rnorm) then
!
!           %---------------------------------------%
!           | No need for further refinement.       |
!           | The cosine of the angle between the   |
!           | corrected residual vector and the old |
!           | residual vector is greater than 0.717 |
!           | In other words the corrected residual |
!           | and the old residual vector share an  |
!           | angle of less than arcCOS(0.717)      |
!           %---------------------------------------%
!
            rnorm = rnorm1
!
         else
!
!           %-------------------------------------------%
!           | Another step of iterative refinement step |
!           | is required. NITREF is used by stat.h     |
!           %-------------------------------------------%
!
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter <= 1) go to 80
!
!           %-------------------------------------------------%
!           | Otherwise RESID is numerically in the span of V |
!           %-------------------------------------------------%
!
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
!
!        %----------------------------------------------%
!        | Branch here directly if iterative refinement |
!        | wasn't necessary or after at most NITER_REF  |
!        | steps of iterative refinement.               |
!        %----------------------------------------------%
!
  100    continue
!
         rstart = .false.
         orth2  = .false.
!
         call cpu_time (t5)
         titref = titref + (t5 - t4)
!
!        %------------------------------------%
!        | STEP 6: Update  j = j+1;  Continue |
!        %------------------------------------%
!
         j = j + 1
         if (j > k+np) then
            call cpu_time (t1)
            tnaitr = tnaitr + (t1 - t0)
            ido = 99
            do 110 i = max(1,k), k+np-1
!
!              %--------------------------------------------%
!              | Check for splitting and deflation.         |
!              | Use a standard test as in the QR algorithm |
!              | REFERENCE: LAPACK subroutine slahqr        |
!              %--------------------------------------------%
!
               tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
               if( tst1==zero ) &
                    tst1 = slanhs( '1', k+np, h, ldh, workd(n+1) )
               if( abs( h( i+1,i ) )<=max( ulp*tst1, smlnum ) ) &
                    h(i+1,i) = zero
 110        continue
!
            if (msglvl > 2) then
               call smout (logfil, k+np, k+np, h, ldh, ndigit, &
                '_naitr: Final upper Hessenberg matrix H of order K+NP')
            end if
!
            go to 9000
         end if
!
!        %--------------------------------------------------------%
!        | Loop back to extend the factorization by another step. |
!        %--------------------------------------------------------%
!
      go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 9000 continue
      return
!
!     %---------------%
!     | End of snaitr |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: snapps
!
!\Description:
!  Given the Arnoldi factorization
!
!     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!
!  apply NP implicit shifts resulting in
!
!     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!
!  where Q is an orthogonal matrix which is the product of rotations
!  and reflections resulting from the NP bulge chage sweeps.
!  The updated Arnoldi factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call snapps
!     ( N, KEV, NP, SHIFTR, SHIFTI, V, LDV, H, LDH, RESID, Q, LDQ,
!       WORKL, WORKD )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Problem size, i.e. size of matrix A.
!
!  KEV     Integer.  (INPUT/OUTPUT)
!          KEV+NP is the size of the input matrix H.
!          KEV is the size of the updated matrix HNEW.  KEV is only
!          updated on ouput when fewer than NP shifts are applied in
!          order to keep the conjugate pair together.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFTR, Real array of length NP.  (INPUT)
!  SHIFTI  Real and imaginary part of the shifts to be applied.
!          Upon, entry to snapps, the shifts must be sorted so that the
!          conjugate pairs are in consecutive locations.
!
!  V       Real N by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, V contains the current KEV+NP Arnoldi vectors.
!          On OUTPUT, V contains the updated KEV Arnoldi vectors
!          in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
!          On INPUT, H contains the current KEV+NP by KEV+NP upper
!          Hessenber matrix of the Arnoldi factorization.
!          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
!          matrix in the KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          On INPUT, RESID contains the the residual vector r_{k+p}.
!          On OUTPUT, RESID is the update residual vector rnew_{k}
!          in the first KEV locations.
!
!  Q       Real KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations and reflections
!          during the bulge chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Real work array of length (KEV+NP).  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  WORKD   Real work array of length 2*N.  (WORKSPACE)
!          Distributed array used in the application of the accumulated
!          orthogonal matrix Q.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers.
!     smout   ARPACK utility routine that prints matrices.
!     svout   ARPACK utility routine that prints vectors.
!     slabad  LAPACK routine that computes machine constants.
!     slacpy  LAPACK matrix copy routine.
!     slamch  LAPACK routine that determines machine constants.
!     slanhs  LAPACK routine that computes various norms of a matrix.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     slarf   LAPACK routine that applies Householder reflection to
!             a matrix.
!     slarfg  LAPACK Householder reflection construction routine.
!     slartg  LAPACK Givens rotation construction routine.
!     slaset  LAPACK matrix initialization routine.
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     saxpy   Level 1 BLAS that computes a vector triad.
!     scopy   Level 1 BLAS that copies one vector to another .
!     sscal   Level 1 BLAS that scales a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: napps.F   SID: 2.4   DATE OF SID: 3/28/97   RELEASE: 2
!
!\Remarks
!  1. In this version, each shift is applied to all the sublocks of
!     the Hessenberg matrix H and not just to the submatrix that it
!     comes from. Deflation as in LAPACK routine slahqr (QR algorithm
!     for upper Hessenberg matrices ) is used.
!     The subdiagonals of H are enforced to be non-negative.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine snapps &
         ( n, kev, np, shiftr, shifti, v, ldv, h, ldh, resid, q, ldq, &
           workl, workd )
!
!! SNAPPS applies implicit shifts to the Arnoldi factorization.
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    kev, ldh, ldq, ldv, n, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 h(ldh,kev+np), resid(n), shifti(np), shiftr(np), &
                 v(ldv,kev+np), q(ldq,kev+np), workd(2*n), workl(kev+np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      integer    i, iend, ir, istart, j, jj, kplusp, msglvl, nr
      logical    cconj, first
      Real &
                 c, f, g, h11, h12, h21, h22, h32, ovfl, r, s, sigmai, &
                 sigmar, smlnum, ulp, unfl, u(3), t, tau, tst1
      save       first, ovfl, smlnum, ulp, unfl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   saxpy, scopy, sscal, slacpy, slarfg, slarf, &
                 slaset, slabad, slartg
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slamch, slanhs, slapy2
      external   slamch, slanhs, slapy2
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs, max, min
!
!     %----------------%
!     | Data statments |
!     %----------------%
!
      data       first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (first) then
!
!        %-----------------------------------------------%
!        | Set machine-dependent constants for the       |
!        | stopping criterion. If norm(H) <= sqrt(OVFL), |
!        | overflow should not occur.                    |
!        | REFERENCE: LAPACK subroutine slahqr           |
!        %-----------------------------------------------%
!
         unfl = slamch( 'safe minimum' )
         ovfl = one / unfl
         call slabad( unfl, ovfl )
         ulp = slamch( 'precision' )
         smlnum = unfl*( n / ulp )
         first = .false.
      end if
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = mnapps
      kplusp = kev + np
!
!     %--------------------------------------------%
!     | Initialize Q to the identity to accumulate |
!     | the rotations and reflections              |
!     %--------------------------------------------%
!
      call slaset ('All', kplusp, kplusp, zero, one, q, ldq)
!
!     %----------------------------------------------%
!     | Quick return if there are no shifts to apply |
!     %----------------------------------------------%
!
      if (np == 0) go to 9000
!
!     %----------------------------------------------%
!     | Chase the bulge with the application of each |
!     | implicit shift. Each shift is applied to the |
!     | whole matrix including each block.           |
!     %----------------------------------------------%
!
      cconj = .false.
      do 110 jj = 1, np
         sigmar = shiftr(jj)
         sigmai = shifti(jj)
!
         if (msglvl > 2 ) then
            call ivout (logfil, 1, jj, ndigit, &
                     '_napps: shift number.')
            call svout (logfil, 1, sigmar, ndigit, &
                     '_napps: The real part of the shift ')
            call svout (logfil, 1, sigmai, ndigit, &
                     '_napps: The imaginary part of the shift ')
         end if
!
!        %-------------------------------------------------%
!        | The following set of conditionals is necessary  |
!        | in order that complex conjugate pairs of shifts |
!        | are applied together or not at all.             |
!        %-------------------------------------------------%
!
         if ( cconj ) then
!
!           %-----------------------------------------%
!           | cconj = .true. means the previous shift |
!           | had non-zero imaginary part.            |
!           %-----------------------------------------%
!
            cconj = .false.
            go to 110
         else if ( jj < np .and. abs( sigmai ) > zero ) then
!
!           %------------------------------------%
!           | Start of a complex conjugate pair. |
!           %------------------------------------%
!
            cconj = .true.
         else if ( jj == np .and. abs( sigmai ) > zero ) then
!
!           %----------------------------------------------%
!           | The last shift has a nonzero imaginary part. |
!           | Don't apply it; thus the order of the        |
!           | compressed H is order KEV+1 since only np-1  |
!           | were applied.                                |
!           %----------------------------------------------%
!
            kev = kev + 1
            go to 110
         end if
         istart = 1
   20    continue
!
!        %--------------------------------------------------%
!        | if sigmai = 0 then                               |
!        |    Apply the jj-th shift ...                     |
!        | else                                             |
!        |    Apply the jj-th and (jj+1)-th together ...    |
!        |    (Note that jj < np at this point in the code) |
!        | end                                              |
!        | to the current block of H. The next do loop      |
!        | determines the current block ;                   |
!        %--------------------------------------------------%
!
         do 30 i = istart, kplusp-1
!
!           %----------------------------------------%
!           | Check for splitting and deflation. Use |
!           | a standard test as in the QR algorithm |
!           | REFERENCE: LAPACK subroutine slahqr    |
!           %----------------------------------------%
!
            tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
            if( tst1==zero ) &
               tst1 = slanhs( '1', kplusp-jj+1, h, ldh, workl )
            if( abs( h( i+1,i ) )<=max( ulp*tst1, smlnum ) ) then
               if (msglvl > 0) then
                  call ivout (logfil, 1, i, ndigit, &
                       '_napps: matrix splitting at row/column no.')
                  call ivout (logfil, 1, jj, ndigit, &
                       '_napps: matrix splitting with shift number.')
                  call svout (logfil, 1, h(i+1,i), ndigit, &
                       '_napps: off diagonal element.')
               end if
               iend = i
               h(i+1,i) = zero
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
!
         if (msglvl > 2) then
             call ivout (logfil, 1, istart, ndigit, &
                         '_napps: Start of current block ')
             call ivout (logfil, 1, iend, ndigit, &
                         '_napps: End of current block ')
         end if
!
!        %------------------------------------------------%
!        | No reason to apply a shift to block of order 1 |
!        %------------------------------------------------%
!
         if ( istart == iend ) go to 100
!
!        %------------------------------------------------------%
!        | If istart + 1 = iend then no reason to apply a       |
!        | complex conjugate pair of shifts on a 2 by 2 matrix. |
!        %------------------------------------------------------%
!
         if ( istart + 1 == iend .and. abs( sigmai ) > zero ) &
            go to 100
!
         h11 = h(istart,istart)
         h21 = h(istart+1,istart)
         if ( abs( sigmai ) <= zero ) then
!
!           %---------------------------------------------%
!           | Real-valued shift ==> apply single shift QR |
!           %---------------------------------------------%
!
            f = h11 - sigmar
            g = h21
!
            do 80 i = istart, iend-1
!
!              %-----------------------------------------------------%
!              | Contruct the plane rotation G to zero out the bulge |
!              %-----------------------------------------------------%
!
               call slartg (f, g, c, s, r)
               if (i > istart) then
!
!                 %-------------------------------------------%
!                 | The following ensures that h(1:iend-1,1), |
!                 | the first iend-2 off diagonal of elements |
!                 | H, remain non negative.                   |
!                 %-------------------------------------------%
!
                  if (r < zero) then
                     r = -r
                     c = -c
                     s = -s
                  end if
                  h(i,i-1) = r
                  h(i+1,i-1) = zero
               end if
!
!              %---------------------------------------------%
!              | Apply rotation to the left of H;  H <- G'*H |
!              %---------------------------------------------%
!
               do 50 j = i, kplusp
                  t        =  c*h(i,j) + s*h(i+1,j)
                  h(i+1,j) = -s*h(i,j) + c*h(i+1,j)
                  h(i,j)   = t
   50          continue
!
!              %---------------------------------------------%
!              | Apply rotation to the right of H;  H <- H*G |
!              %---------------------------------------------%
!
               do j = 1, min(i+2,iend)
                  t        =  c*h(j,i) + s*h(j,i+1)
                  h(j,i+1) = -s*h(j,i) + c*h(j,i+1)
                  h(j,i)   = t
               end do
!
!              %----------------------------------------------------%
!              | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!              %----------------------------------------------------%
!
               do j = 1, min( i+jj, kplusp )
                  t        =   c*q(j,i) + s*q(j,i+1)
                  q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                  q(j,i)   = t
               end do
!
!              %---------------------------%
!              | Prepare for next rotation |
!              %---------------------------%
!
               if (i < iend-1) then
                  f = h(i+1,i)
                  g = h(i+2,i)
               end if
   80       continue
!
!           %-----------------------------------%
!           | Finished applying the real shift. |
!           %-----------------------------------%
!
         else
!
!           %----------------------------------------------------%
!           | Complex conjugate shifts ==> apply double shift QR |
!           %----------------------------------------------------%
!
            h12 = h(istart,istart+1)
            h22 = h(istart+1,istart+1)
            h32 = h(istart+2,istart+1)
!
!           %---------------------------------------------------------%
!           | Compute 1st column of (H - shift*I)*(H - conj(shift)*I) |
!           %---------------------------------------------------------%
!
            s    = 2.0*sigmar
            t = slapy2 ( sigmar, sigmai )
            u(1) = ( h11 * (h11 - s) + t * t ) / h21 + h12
            u(2) = h11 + h22 - s
            u(3) = h32
!
            do 90 i = istart, iend-1
!
               nr = min ( 3, iend-i+1 )
!
!              %-----------------------------------------------------%
!              | Construct Householder reflector G to zero out u(1). |
!              | G is of the form I - tau*( 1 u )' * ( 1 u' ).       |
!              %-----------------------------------------------------%
!
               call slarfg ( nr, u(1), u(2), 1, tau )
!
               if (i > istart) then
                  h(i,i-1)   = u(1)
                  h(i+1,i-1) = zero
                  if (i < iend-1) h(i+2,i-1) = zero
               end if
               u(1) = one
!
!              %--------------------------------------%
!              | Apply the reflector to the left of H |
!              %--------------------------------------%
!
               call slarf ('Left', nr, kplusp-i+1, u, 1, tau, &
                           h(i,i), ldh, workl)
!
!              %---------------------------------------%
!              | Apply the reflector to the right of H |
!              %---------------------------------------%
!
               ir = min ( i+3, iend )
               call slarf ('Right', ir, nr, u, 1, tau, &
                           h(1,i), ldh, workl)
!
!              %-----------------------------------------------------%
!              | Accumulate the reflector in the matrix Q;  Q <- Q*G |
!              %-----------------------------------------------------%
!
               call slarf ('Right', kplusp, nr, u, 1, tau, &
                           q(1,i), ldq, workl)
!
!              %----------------------------%
!              | Prepare for next reflector |
!              %----------------------------%
!
               if (i < iend-1) then
                  u(1) = h(i+1,i)
                  u(2) = h(i+2,i)
                  if (i < iend-2) u(3) = h(i+3,i)
               end if
!
   90       continue
!
!           %--------------------------------------------%
!           | Finished applying a complex pair of shifts |
!           | to the current block                       |
!           %--------------------------------------------%
!
         end if
!
  100    continue
!
!        %---------------------------------------------------------%
!        | Apply the same shift to the next block if there is any. |
!        %---------------------------------------------------------%
!
         istart = iend + 1
         if (iend < kplusp) go to 20
!
!        %---------------------------------------------%
!        | Loop back to the top to get the next shift. |
!        %---------------------------------------------%
!
  110 continue
!
!     %--------------------------------------------------%
!     | Perform a similarity transformation that makes   |
!     | sure that H will have non negative sub diagonals |
!     %--------------------------------------------------%
!
      do 120 j=1,kev
         if ( h(j+1,j) < zero ) then
              call sscal( kplusp-j+1, -one, h(j+1,j), ldh )
              call sscal( min(j+2, kplusp), -one, h(1,j+1), 1 )
              call sscal( min(j+np+1,kplusp), -one, q(1,j+1), 1 )
         end if
 120  continue
!
      do 130 i = 1, kev
!
!        %--------------------------------------------%
!        | Final check for splitting and deflation.   |
!        | Use a standard test as in the QR algorithm |
!        | REFERENCE: LAPACK subroutine slahqr        |
!        %--------------------------------------------%
!
         tst1 = abs( h( i, i ) ) + abs( h( i+1, i+1 ) )
         if( tst1==zero ) &
             tst1 = slanhs( '1', kev, h, ldh, workl )
         if( h( i+1,i ) <= max( ulp*tst1, smlnum ) ) &
             h(i+1,i) = zero
 130  continue
!
!     %-------------------------------------------------%
!     | Compute the (kev+1)-st column of (V*Q) and      |
!     | temporarily store the result in WORKD(N+1:2*N). |
!     | This is needed in the residual update since we  |
!     | cannot GUARANTEE that the corresponding entry   |
!     | of H would be zero as in exact arithmetic.      |
!     %-------------------------------------------------%
!
      if (h(kev+1,kev) > zero) &
          call sgemv ('N', n, kplusp, one, v, ldv, q(1,kev+1), 1, zero, &
                      workd(n+1), 1)
!
!     %----------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order       |
!     | taking advantage of the upper Hessenberg structure of Q. |
!     %----------------------------------------------------------%
!
      do 140 i = 1, kev
         call sgemv ('N', n, kplusp-i+1, one, v, ldv, &
                     q(1,kev-i+1), 1, zero, workd, 1)
         call scopy (n, workd, 1, v(1,kplusp-i+1), 1)
  140 continue
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
      call slacpy ('A', n, kev, v(1,kplusp-kev+1), ldv, v, ldv)
!
!     %--------------------------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
!     %--------------------------------------------------------------%
!
      if (h(kev+1,kev) > zero) &
         call scopy (n, workd(n+1), 1, v(1,kev+1), 1)
!
!     %-------------------------------------%
!     | Update the residual vector:         |
!     |    r <- sigmak*r + betak*v(:,kev+1) |
!     | where                               |
!     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
!     |    betak = e_{kev+1}'*H*e_{kev}     |
!     %-------------------------------------%
!
      call sscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,kev) > zero) &
         call saxpy (n, h(kev+1,kev), v(1,kev+1), 1, resid, 1)
!
      if (msglvl > 1) then
         call svout (logfil, 1, q(kplusp,kev), ndigit, &
              '_napps: sigmak = (e_{kev+p}^T*Q)*e_{kev}')
         call svout (logfil, 1, h(kev+1,kev), ndigit, &
              '_napps: betak = e_{kev+1}^T*H*e_{kev}')
         call ivout (logfil, 1, kev, ndigit, &
                     '_napps: Order of the final Hessenberg matrix ')
         if (msglvl > 2) then
            call smout (logfil, kev, kev, h, ldh, ndigit, &
            '_napps: updated Hessenberg matrix H for next iteration')
         end if
!
      end if
!
 9000 continue
      call cpu_time (t1)
      tnapps = tnapps + (t1 - t0)
!
      return
!
!     %---------------%
!     | End of snapps |
!     %---------------%
!
      end
!\BeginDoc
!
!\Name: snaup2
!
!\Description:
!  Intermediate level interface called by snaupd.
!
!\Usage:
!  call snaup2
!     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZR, RITZI, BOUNDS,
!       Q, LDQ, WORKL, IPNTR, WORKD, INFO )
!
!\Arguments
!
!  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in snaupd.
!  MODE, ISHIFT, MXITER: see the definition of IPARAM in snaupd.
!
!  NP      Integer.  (INPUT/OUTPUT)
!          Contains the number of implicit shifts to apply during
!          each Arnoldi iteration.
!          If ISHIFT=1, NP is adjusted dynamically at each iteration
!          to accelerate convergence and prevent stagnation.
!          This is also roughly equal to the number of matrix-vector
!          products (involving the operator OP) per Arnoldi iteration.
!          The logic for adjusting is contained within the current
!          subroutine.
!          If ISHIFT=0, NP is the number of shifts the user needs
!          to provide via reverse comunication. 0 < NP < NCV-NEV.
!          NP may be less than NCV-NEV for two reasons. The first, is
!          to keep complex conjugate pairs of "wanted" Ritz values
!          together. The second, is that a leading block of the current
!          upper Hessenberg matrix has split off and contains "unwanted"
!          Ritz values.
!          Upon termination of the IRA iteration, NP contains the number
!          of "converged" wanted Ritz values.
!
!  IUPD    Integer.  (INPUT)
!          IUPD .EQ. 0: use explicit restart instead implicit update.
!          IUPD .NE. 0: use implicit update.
!
!  V       Real N by (NEV+NP) array.  (INPUT/OUTPUT)
!          The Arnoldi basis vectors are returned in the first NEV
!          columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (NEV+NP) by (NEV+NP) array.  (OUTPUT)
!          H is used to store the generated upper Hessenberg matrix
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZR,  Real arrays of length NEV+NP.  (OUTPUT)
!  RITZI   RITZR(1:NEV) (resp. RITZI(1:NEV)) contains the real (resp.
!          imaginary) part of the computed Ritz values of OP.
!
!  BOUNDS  Real array of length NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) contain the error bounds corresponding to
!          the computed Ritz values.
!
!  Q       Real (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
!          Private (replicated) work array used to accumulate the
!          rotation in the shift application step.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Real work array of length at least
!          (NEV+NP)**2 + 3*(NEV+NP).  (INPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  It is used in shifts calculation, shifts
!          application and convergence checking.
!
!          On exit, the last 3*(NEV+NP) locations of WORKL contain
!          the Ritz values (real,imaginary) and associated Ritz
!          estimates of the current Hessenberg matrix.  They are
!          listed in the same order as returned from sneigh.
!
!          If ISHIFT .EQ. O and IDO .EQ. 3, the first 2*NP locations
!          of WORKL are used in reverse communication to hold the user
!          supplied shifts.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD for
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!
!  WORKD   Real work array of length 3*N.  (WORKSPACE)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration !!!!!!!!!!
!          See Data Distribution Note in DNAUPD.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =     0: Normal return.
!          =     1: Maximum number of iterations taken.
!                   All possible eigenvalues of OP has been found.
!                   NP returns the number of converged Ritz values.
!          =     2: No shifts could be applied.
!          =    -8: Error return from LAPACK eigenvalue calculation;
!                   This should never happen.
!          =    -9: Starting vector is zero.
!          = -9999: Could not build an Arnoldi factorization.
!                   Size that was built in returned in NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     sgetv0  ARPACK initial vector generation routine.
!     snaitr  ARPACK Arnoldi factorization routine.
!     snapps  ARPACK application of implicit shifts routine.
!     snconv  ARPACK convergence of Ritz values routine.
!     sneigh  ARPACK compute Ritz values and error bounds routine.
!     sngets  ARPACK reorder Ritz values and error bounds routine.
!     ssortc  ARPACK sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     smout   ARPACK utility routine that prints matrices
!     svout   ARPACK utility routine that prints vectors.
!     slamch  LAPACK routine that determines machine constants.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     scopy   Level 1 BLAS that copies one vector to another .
!     sdot    Level 1 BLAS that computes the scalar product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sswap   Level 1 BLAS that swaps two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: naup2.F   SID: 2.8   DATE OF SID: 10/17/00   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine snaup2 &
         ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, &
           ishift, mxiter, v, ldv, h, ldh, ritzr, ritzi, bounds, &
           q, ldq, workl, ipntr, workd, info )
!
!! SNAUP2 is an intermediate interface called by SNAUPD.
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, mode, ldh, ldq, ldv, mxiter, &
                 n, nev, np
      Real &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(13)
      Real &
                 bounds(nev+np), h(ldh,nev+np), q(ldq,nev+np), resid(n), &
                 ritzi(nev+np), ritzr(nev+np), v(ldv,nev+np), &
                 workd(3*n), workl( (nev+np)*(nev+np+3) )
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  wprime*2
      logical    cnorm , getv0, initv, update, ushift
      integer    ierr  , iter , j    , kplusp, msglvl, nconv, &
                 nevbef, nev0 , np0  , nptemp, numcnv
      Real &
                 rnorm , temp , eps23
      save       cnorm , getv0, initv, update, ushift, &
                 rnorm , iter , eps23, kplusp, msglvl, nconv , &
                 nevbef, nev0 , np0  , numcnv
!
!     %-----------------------%
!     | Local array arguments |
!     %-----------------------%
!
      integer    kp(4)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy , sgetv0, snaitr, snconv, sneigh, &
                 sngets, snapps, svout , ivout 
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 sdot, snrm2, slapy2, slamch
      external   sdot, snrm2, slapy2, slamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min, max, abs, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido == 0) then
!
         call cpu_time (t0)
!
         msglvl = mnaup2
!
!        %-------------------------------------%
!        | Get the machine dependent constant. |
!        %-------------------------------------%
!
         eps23 = slamch('Epsilon-Machine')
         eps23 = eps23**(2.0E+0 / 3.0E+0)
!
         nev0   = nev
         np0    = np
!
!        %-------------------------------------%
!        | kplusp is the bound on the largest  |
!        |        Lanczos factorization built. |
!        | nconv is the current number of      |
!        |        "converged" eigenvlues.      |
!        | iter is the counter on the current  |
!        |      iteration step.                |
!        %-------------------------------------%
!
         kplusp = nev + np
         nconv  = 0
         iter   = 0
!
!        %---------------------------------------%
!        | Set flags for computing the first NEV |
!        | steps of the Arnoldi factorization.   |
!        %---------------------------------------%
!
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
!
         if (info /= 0) then
!
!           %--------------------------------------------%
!           | User provides the initial residual vector. |
!           %--------------------------------------------%
!
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
!
!     %---------------------------------------------%
!     | Get a possibly random starting vector and   |
!     | force it into the range of the operator OP. |
!     %---------------------------------------------%
!
   10 continue
!
      if (getv0) then
         call sgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm, &
                      ipntr, workd, info)
!
         if (ido /= 99) go to 9000
!
         if (rnorm == zero) then
!
!           %-----------------------------------------%
!           | The initial vector is zero. Error exit. |
!           %-----------------------------------------%
!
            info = -9
            go to 1100
         end if
         getv0 = .false.
         ido  = 0
      end if
!
!     %-----------------------------------%
!     | Back from reverse communication : |
!     | continue with update step         |
!     %-----------------------------------%
!
      if (update) go to 20
!
!     %-------------------------------------------%
!     | Back from computing user specified shifts |
!     %-------------------------------------------%
!
      if (ushift) go to 50
!
!     %-------------------------------------%
!     | Back from computing residual norm   |
!     | at the end of the current iteration |
!     %-------------------------------------%
!
      if (cnorm)  go to 100
!
!     %----------------------------------------------------------%
!     | Compute the first NEV steps of the Arnoldi factorization |
!     %----------------------------------------------------------%
!
      call snaitr (ido, bmat, n, 0, nev, mode, resid, rnorm, v, ldv, &
                   h, ldh, ipntr, workd, info)
!
!     %---------------------------------------------------%
!     | ido /= 99 implies use of reverse communication  |
!     | to compute operations involving OP and possibly B |
!     %---------------------------------------------------%
!
      if (ido /= 99) go to 9000

      if (info > 0) then
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
!     |           Each iteration implicitly restarts the Arnoldi     |
!     |           factorization in place.                            |
!     |                                                              |
!     %--------------------------------------------------------------%
!
 1000 continue
!
         iter = iter + 1
!
         if (msglvl > 0) then
            call ivout (logfil, 1, iter, ndigit, &
                 '_naup2: **** Start of major iteration number ****')
         end if
!
!        %-----------------------------------------------------------%
!        | Compute NP additional steps of the Arnoldi factorization. |
!        | Adjust NP since NEV might have been updated by last call  |
!        | to the shift application routine snapps.                  |
!        %-----------------------------------------------------------%
!
         np  = kplusp - nev
!
         if (msglvl > 1) then
            call ivout (logfil, 1, nev, ndigit, &
           '_naup2: The length of the current Arnoldi factorization')
            call ivout (logfil, 1, np, ndigit, &
                 '_naup2: Extend the Arnoldi factorization by')
         end if
!
!        %-----------------------------------------------------------%
!        | Compute NP additional steps of the Arnoldi factorization. |
!        %-----------------------------------------------------------%
!
         ido = 0
   20    continue
         update = .true.

         call snaitr (ido  , bmat, n  , nev, np , mode , resid, &
                      rnorm, v   , ldv, h  , ldh, ipntr, workd, &
                      info)
!
!        %---------------------------------------------------%
!        | ido /= 99 implies use of reverse communication  |
!        | to compute operations involving OP and possibly B |
!        %---------------------------------------------------%
!
         if (ido /= 99) go to 9000

         if (info > 0) then
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.

         if (msglvl > 1) then
            call svout (logfil, 1, rnorm, ndigit, &
                 '_naup2: Corresponding B-norm of the residual')
         end if
!
!        %--------------------------------------------------------%
!        | Compute the eigenvalues and corresponding error bounds |
!        | of the current upper Hessenberg matrix.                |
!        %--------------------------------------------------------%
!
         call sneigh (rnorm, kplusp, h, ldh, ritzr, ritzi, bounds, &
                      q, ldq, workl, ierr)
!
         if (ierr /= 0) then
            info = -8
            go to 1200
         end if
!
!        %----------------------------------------------------%
!        | Make a copy of eigenvalues and corresponding error |
!        | bounds obtained from sneigh.                       |
!        %----------------------------------------------------%
!
         call scopy(kplusp, ritzr, 1, workl(kplusp**2+1), 1)
         call scopy(kplusp, ritzi, 1, workl(kplusp**2+kplusp+1), 1)
         call scopy(kplusp, bounds, 1, workl(kplusp**2+2*kplusp+1), 1)
!
!        %---------------------------------------------------%
!        | Select the wanted Ritz values and their bounds    |
!        | to be used in the convergence test.               |
!        | The wanted part of the spectrum and corresponding |
!        | error bounds are in the last NEV loc. of RITZR,   |
!        | RITZI and BOUNDS respectively. The variables NEV  |
!        | and NP may be updated if the NEV-th wanted Ritz   |
!        | value has a non zero imaginary part. In this case |
!        | NEV is increased by one and NP decreased by one.  |
!        | NOTE: The last two arguments of sngets are no     |
!        | longer used as of version 2.1.                    |
!        %---------------------------------------------------%
!
         nev = nev0
         np = np0
         numcnv = nev
         call sngets (ishift, which, nev, np, ritzr, ritzi, &
                      bounds, workl, workl(np+1))
         if (nev == nev0+1) numcnv = nev0+1
!
!        %-------------------%
!        | Convergence test. |
!        %-------------------%
!
         call scopy (nev, bounds(np+1), 1, workl(2*np+1), 1)
         call snconv (nev, ritzr(np+1), ritzi(np+1), workl(2*np+1), &
              tol, nconv)
!
         if (msglvl > 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = numcnv
            kp(4) = nconv
            call ivout (logfil, 4, kp, ndigit, &
                        '_naup2: NEV, NP, NUMCNV, NCONV are')
            call svout (logfil, kplusp, ritzr, ndigit, &
                 '_naup2: Real part of the eigenvalues of H')
            call svout (logfil, kplusp, ritzi, ndigit, &
                 '_naup2: Imaginary part of the eigenvalues of H')
            call svout (logfil, kplusp, bounds, ndigit, &
                '_naup2: Ritz estimates of the current NCV Ritz values')
         end if
!
!        %---------------------------------------------------------%
!        | Count the number of unwanted Ritz values that have zero |
!        | Ritz estimates. If any Ritz estimates are equal to zero |
!        | then a leading block of H of order equal to at least    |
!        | the number of Ritz values with zero Ritz estimates has  |
!        | split off. None of these Ritz values may be removed by  |
!        | shifting. Decrease NP the number of shifts to apply. If |
!        | no shifts may be applied, then prepare to exit          |
!        %---------------------------------------------------------%
!
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) == zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
!
         if ( (nconv >= numcnv) .or. &
              (iter > mxiter) .or. &
              (np == 0) ) then
!
            if (msglvl > 4) then
               call svout(logfil, kplusp, workl(kplusp**2+1), ndigit, &
                   '_naup2: Real part of the eig computed by _neigh:')
               call svout(logfil, kplusp, workl(kplusp**2+kplusp+1), &
                           ndigit, &
                   '_naup2: Imag part of the eig computed by _neigh:')
               call svout(logfil, kplusp, workl(kplusp**2+kplusp*2+1), &
                           ndigit, &
                   '_naup2: Ritz eistmates computed by _neigh:')
            end if
!
!           %------------------------------------------------%
!           | Prepare to exit. Put the converged Ritz values |
!           | and corresponding bounds in RITZ(1:NCONV) and  |
!           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
!           | careful when NCONV > NP                        |
!           %------------------------------------------------%
!
!           %------------------------------------------%
!           |  Use h( 3,1 ) as storage to communicate  |
!           |  rnorm to _neupd if needed               |
!           %------------------------------------------%

            h(3,1) = rnorm
!
!           %----------------------------------------------%
!           | To be consistent with sngets, we first do a  |
!           | pre-processing sort in order to keep complex |
!           | conjugate pairs together.  This is similar   |
!           | to the pre-processing sort used in sngets    |
!           | except that the sort is done in the opposite |
!           | order.                                       |
!           %----------------------------------------------%
!
            if (which == 'LM') wprime = 'SR'
            if (which == 'SM') wprime = 'LR'
            if (which == 'LR') wprime = 'SM'
            if (which == 'SR') wprime = 'LM'
            if (which == 'LI') wprime = 'SM'
            if (which == 'SI') wprime = 'LM'
!
            call ssortc (wprime, .true., kplusp, ritzr, ritzi, bounds)
!
!           %----------------------------------------------%
!           | Now sort Ritz values so that converged Ritz  |
!           | values appear within the first NEV locations |
!           | of ritzr, ritzi and bounds, and the most     |
!           | desired one appears at the front.            |
!           %----------------------------------------------%
!
            if (which == 'LM') wprime = 'SM'
            if (which == 'SM') wprime = 'LM'
            if (which == 'LR') wprime = 'SR'
            if (which == 'SR') wprime = 'LR'
            if (which == 'LI') wprime = 'SI'
            if (which == 'SI') wprime = 'LI'
!
            call ssortc(wprime, .true., kplusp, ritzr, ritzi, bounds)
!
!           %--------------------------------------------------%
!           | Scale the Ritz estimate of each Ritz value       |
!           | by 1 / max(eps23,magnitude of the Ritz value).   |
!           %--------------------------------------------------%
!
            do 35 j = 1, numcnv
                temp = max(eps23,slapy2(ritzr(j), &
                                         ritzi(j)))
                bounds(j) = bounds(j)/temp
 35         continue
!
!           %----------------------------------------------------%
!           | Sort the Ritz values according to the scaled Ritz  |
!           | esitmates.  This will push all the converged ones  |
!           | towards the front of ritzr, ritzi, bounds          |
!           | (in the case when NCONV < NEV.)                    |
!           %----------------------------------------------------%
!
            wprime = 'LR'
            call ssortc(wprime, .true., numcnv, bounds, ritzr, ritzi)
!
!           %----------------------------------------------%
!           | Scale the Ritz estimate back to its original |
!           | value.                                       |
!           %----------------------------------------------%
!
            do j = 1, numcnv
                temp = max(eps23, slapy2(ritzr(j), &
                                         ritzi(j)))
                bounds(j) = bounds(j)*temp
            end do
!
!           %------------------------------------------------%
!           | Sort the converged Ritz values again so that   |
!           | the "threshold" value appears at the front of  |
!           | ritzr, ritzi and bound.                        |
!           %------------------------------------------------%
!
            call ssortc(which, .true., nconv, ritzr, ritzi, bounds)
!
            if (msglvl > 1) then
               call svout (logfil, kplusp, ritzr, ndigit, &
                  '_naup2: Sorted real part of the eigenvalues')
               call svout (logfil, kplusp, ritzi, ndigit, &
                  '_naup2: Sorted imaginary part of the eigenvalues')
               call svout (logfil, kplusp, bounds, ndigit, &
                  '_naup2: Sorted ritz estimates.')
            end if
!
!           %------------------------------------%
!           | Max iterations have been exceeded. |
!           %------------------------------------%
!
            if (iter > mxiter .and. nconv < numcnv) info = 1
!
!           %---------------------%
!           | No shifts to apply. |
!           %---------------------%
!
            if (np == 0 .and. nconv < numcnv) info = 2
!
            np = nconv
            go to 1100
!
         else if ( (nconv < numcnv) .and. (ishift == 1) ) then
!
!           %-------------------------------------------------%
!           | Do not have all the requested eigenvalues yet.  |
!           | To prevent possible stagnation, adjust the size |
!           | of NEV.                                         |
!           %-------------------------------------------------%
!
            nevbef = nev
            nev = nev + min(nconv, np/2)
            if (nev == 1 .and. kplusp >= 6) then
               nev = kplusp / 2
            else if (nev == 1 .and. kplusp > 3) then
               nev = 2
            end if
            np = kplusp - nev
!
!           %---------------------------------------%
!           | If the size of NEV was just increased |
!           | resort the eigenvalues.               |
!           %---------------------------------------%
!
            if (nevbef < nev) &
               call sngets (ishift, which, nev, np, ritzr, ritzi, &
                    bounds, workl, workl(np+1))
!
         end if
!
         if (msglvl > 0) then
            call ivout (logfil, 1, nconv, ndigit, &
                 '_naup2: no. of "converged" Ritz values at this iter.')
            if (msglvl > 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, &
                    '_naup2: NEV and NP are')
               call svout (logfil, nev, ritzr(np+1), ndigit, &
                    '_naup2: "wanted" Ritz values -- real part')
               call svout (logfil, nev, ritzi(np+1), ndigit, &
                    '_naup2: "wanted" Ritz values -- imag part')
               call svout (logfil, nev, bounds(np+1), ndigit, &
                    '_naup2: Ritz estimates of the "wanted" values ')
            end if
         end if
!
         if (ishift == 0) then
!
!           %-------------------------------------------------------%
!           | User specified shifts: reverse comminucation to       |
!           | compute the shifts. They are returned in the first    |
!           | 2*NP locations of WORKL.                              |
!           %-------------------------------------------------------%
!
            ushift = .true.
            ido = 3
            go to 9000
         end if
!
   50    continue
!
!        %------------------------------------%
!        | Back from reverse communication;   |
!        | User specified shifts are returned |
!        | in WORKL(1:2*NP)                   |
!        %------------------------------------%
!
         ushift = .false.
!
         if ( ishift == 0 ) then
!
!            %----------------------------------%
!            | Move the NP shifts from WORKL to |
!            | RITZR, RITZI to free up WORKL    |
!            | for non-exact shift case.        |
!            %----------------------------------%
!
             call scopy (np, workl,       1, ritzr, 1)
             call scopy (np, workl(np+1), 1, ritzi, 1)
         end if
!
         if (msglvl > 2) then
            call ivout (logfil, 1, np, ndigit, &
                        '_naup2: The number of shifts to apply ')
            call svout (logfil, np, ritzr, ndigit, &
                        '_naup2: Real part of the shifts')
            call svout (logfil, np, ritzi, ndigit, &
                        '_naup2: Imaginary part of the shifts')
            if ( ishift == 1 ) &
                call svout (logfil, np, bounds, ndigit, &
                        '_naup2: Ritz estimates of the shifts')
         end if
!
!        %---------------------------------------------------------%
!        | Apply the NP implicit shifts by QR bulge chasing.       |
!        | Each shift is applied to the whole upper Hessenberg     |
!        | matrix H.                                               |
!        | The first 2*N locations of WORKD are used as workspace. |
!        %---------------------------------------------------------%
!
         call snapps (n, nev, np, ritzr, ritzi, v, ldv, &
                      h, ldh, resid, q, ldq, workl, workd)
!
!        %---------------------------------------------%
!        | Compute the B-norm of the updated residual. |
!        | Keep B*RESID in WORKD(1:N) to be used in    |
!        | the first step of the next call to snaitr.  |
!        %---------------------------------------------%
!
         cnorm = .true.
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*RESID |
!           %----------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
            call scopy (n, resid, 1, workd, 1)
         end if
!
  100    continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(1:N) := B*RESID            |
!        %----------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         if (bmat == 'G') then
            rnorm = sdot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat == 'I') then
            rnorm = snrm2(n, resid, 1)
         end if
         cnorm = .false.
!
         if (msglvl > 2) then
            call svout (logfil, 1, rnorm, ndigit, &
            '_naup2: B-norm of residual for compressed factorization')
            call smout (logfil, nev, nev, h, ldh, ndigit, &
              '_naup2: Compressed upper Hessenberg matrix H')
         end if
!
      go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 1100 continue
!
      mxiter = iter
      nev = numcnv
!
 1200 continue
      ido = 99
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      call cpu_time (t1)
      tnaup2 = t1 - t0
!
 9000 continue
!
!     %---------------%
!     | End of snaup2 |
!     %---------------%
!
      return
      end
      subroutine snaupd &
         ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
           ipntr, workd, workl, lworkl, info )
!
!! SNAUPD is an interface for the Implicitly Restarted Arnoldi iteration.
!
!\BeginDoc
!
!\Name: snaupd
!
!\Description:
!  Reverse communication interface for the Implicitly Restarted Arnoldi
!  iteration. This subroutine computes approximations to a few eigenpairs
!  of a linear operator "OP" with respect to a semi-inner product defined by
!  a symmetric positive semi-definite real matrix B. B may be the identity
!  matrix. NOTE: If the linear operator "OP" is real and symmetric
!  with respect to the real positive semi-definite symmetric matrix B,
!  i.e. B*OP = (OP`)*B, then subroutine ssaupd should be used instead.
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  snaupd is usually called iteratively to solve one of the
!  following problems:
!
!  Mode 1:  A*x = lambda*x.
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*x = lambda*M*x, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)
!
!  Mode 3:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Real_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then
!           amu = 1/2 * [ 1/(lambda-sigma) + 1/(lambda-conjg(sigma)) ].
!           Note: If sigma is real, i.e. imaginary part of sigma is zero;
!                 Real_Part{ inv[A - sigma*M]*M } == inv[A - sigma*M]*M
!                 amu == 1/(lambda-sigma).
!
!  Mode 4:  A*x = lambda*M*x, M symmetric semi-definite
!           ===> OP = Imaginary_Part{ inv[A - sigma*M]*M }  and  B = M.
!           ===> shift-and-invert mode (in real arithmetic)
!           If OP*x = amu*x, then
!           amu = 1/2i * [ 1/(lambda-sigma) - 1/(lambda-conjg(sigma)) ].
!
!  Both mode 3 and 4 give the same enhancement to eigenvalues close to
!  the (complex) shift sigma.  However, as lambda goes to infinity,
!  the operator OP in mode 4 dampens the eigenvalues more strongly than
!  does OP defined in mode 3.
!
!  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
!        should be accomplished either by a direct method
!        using a sparse matrix factorization and solving
!
!           [A - sigma*M]*w = v  or M*w = v,
!
!        or through an iterative method for solving these
!        systems.  If an iterative method is used, the
!        convergence test must be more stringent than
!        the accuracy requirements for the eigenvalue
!        approximations.
!
!\Usage:
!  call snaupd
!     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to snaupd.  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          snaupd with the result.  The operand is given in
!          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    In mode 3 and 4, the vector B * X is already
!                    available in WORKD(ipntr(3)).  It does not
!                    need to be recomputed in forming OP * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO =  3: compute the IPARAM(8) real and imaginary parts
!                    of the shifts where INPTR(14) is the pointer
!                    into WORKL for placing the shifts. See Remark
!                    5 below.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
!          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  WHICH   Character*2.  (INPUT)
!          'LM' -> want the NEV eigenvalues of largest magnitude.
!          'SM' -> want the NEV eigenvalues of smallest magnitude.
!          'LR' -> want the NEV eigenvalues of largest real part.
!          'SR' -> want the NEV eigenvalues of smallest real part.
!          'LI' -> want the NEV eigenvalues of largest imaginary part.
!          'SI' -> want the NEV eigenvalues of smallest imaginary part.
!
!  NEV     Integer.  (INPUT/OUTPUT)
!          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.
!
!  TOL     Real scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
!          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
!          DEFAULT = SLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine SLAMCH).
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          On INPUT:
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V. NCV must satisfy the two
!          inequalities 2 <= NCV-NEV and NCV <= N.
!          This will indicate how many Arnoldi vectors are generated
!          at each iteration.  After the startup phase in which NEV
!          Arnoldi vectors are generated, the algorithm generates
!          approximately NCV-NEV Arnoldi vectors at each subsequent update
!          iteration. Most of the cost in generating each Arnoldi vector is
!          in the matrix-vector operation OP*x.
!          NOTE: 2 <= NCV-NEV in order that complex conjugate pairs of Ritz
!          values are kept together. (See remark 4 below)
!
!  V       Real array N by NCV.  (OUTPUT)
!          Contains the final set of Arnoldi basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are provided by the user via
!                      reverse communication.  The real and imaginary
!                      parts of the NCV eigenvalues of the Hessenberg
!                      matrix H are returned in the part of the WORKL
!                      array corresponding to RITZR and RITZI. See remark
!                      5 below.
!          ISHIFT = 1: exact shifts with respect to the current
!                      Hessenberg matrix H.  This is equivalent to
!                      restarting the iteration with a starting vector
!                      that is a linear combination of approximate Schur
!                      vectors associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = No longer referenced.
!
!          IPARAM(3) = MXITER
!          On INPUT:  maximum number of Arnoldi update iterations allowed.
!          On OUTPUT: actual number of Arnoldi update iterations taken.
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" Ritz values.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used.
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4; See under \Description of snaupd for the
!          four modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), snaupd returns NP, the number
!          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
!          5 below.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.
!
!  IPNTR   Integer array of length 14.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD and WORKL
!          arrays for matrices/vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X in WORKD.
!          IPNTR(2): pointer to the current result vector Y in WORKD.
!          IPNTR(3): pointer to the vector B * X in WORKD when used in
!                    the shift-and-invert mode.
!          IPNTR(4): pointer to the next available location in WORKL
!                    that is untouched by the program.
!          IPNTR(5): pointer to the NCV by NCV upper Hessenberg matrix
!                    H in WORKL.
!          IPNTR(6): pointer to the real part of the ritz value array
!                    RITZR in WORKL.
!          IPNTR(7): pointer to the imaginary part of the ritz value array
!                    RITZI in WORKL.
!          IPNTR(8): pointer to the Ritz estimates in array WORKL associated
!                    with the Ritz values located in RITZR and RITZI in WORKL.
!
!          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.
!
!          Note: IPNTR(9:13) is only referenced by sneupd. See Remark 2 below.
!
!          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
!                     original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     sneupd if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration. Upon termination
!          WORKD(1:N) contains B*RESID(1:N). If an invariant subspace
!          associated with the converged Ritz values is desired, see remark
!          2 below, subroutine sneupd uses this output.
!          See Data Distribution Note below.
!
!  WORKL   Real work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least 3*NCV**2 + 6*NCV.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the
!                Implicitly restarted Arnoldi iteration. One possibility
!                is to increase the size of NCV relative to NEV.
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iteration
!                must be greater than zero.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array is not sufficient.
!          = -8: Error return from LAPACK eigenvalue calculation;
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization.
!
!\Remarks
!  1. The computed Ritz values are approximate eigenvalues of OP. The
!     selection of WHICH should be made with this in mind when
!     Mode = 3 and 4.  After convergence, approximate eigenvalues of the
!     original problem may be obtained with the ARPACK subroutine sneupd.
!
!  2. If a basis for the invariant subspace corresponding to the converged Ritz
!     values is needed, the user must call sneupd immediately following
!     completion of snaupd. This is new starting with release 2 of ARPACK.
!
!  3. If M can be factored into a Cholesky factorization M = LL`
!     then Mode = 2 should not be selected.  Instead one should use
!     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
!     linear systems should be solved with L and L` rather
!     than computing inverses.  After convergence, an approximate
!     eigenvector z of the original problem is recovered by solving
!     L`z = x  where x is a Ritz vector of OP.
!
!  4. At present there is no a-priori analysis to guide the selection
!     of NCV relative to NEV.  The only formal requrement is that NCV > NEV + 2.
!     However, it is recommended that NCV >= 2*NEV+1.  If many problems of
!     the same type are to be solved, one should experiment with increasing
!     NCV while keeping NEV fixed for a given test problem.  This will
!     usually decrease the required number of OP*x operations but it
!     also increases the work and storage required to maintain the orthogonal
!     basis vectors.  The optimal "cross-over" with respect to CPU time
!     is problem dependent and must be determined empirically.
!     See Chapter 8 of Reference 2 for further information.
!
!  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
!     NP = IPARAM(8) real and imaginary parts of the shifts in locations
!         real part                  imaginary part
!         -----------------------    --------------
!     1   WORKL(IPNTR(14))           WORKL(IPNTR(14)+NP)
!     2   WORKL(IPNTR(14)+1)         WORKL(IPNTR(14)+NP+1)
!                        .                          .
!                        .                          .
!                        .                          .
!     NP  WORKL(IPNTR(14)+NP-1)      WORKL(IPNTR(14)+2*NP-1).
!
!     Only complex conjugate pairs of shifts may be applied and the pairs
!     must be placed in consecutive locations. The real part of the
!     eigenvalues of the current upper Hessenberg matrix are located in
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1) and the imaginary part
!     in WORKL(IPNTR(7)) through WORKL(IPNTR(7)+NCV-1). They are ordered
!     according to the order defined by WHICH. The complex conjugate
!     pairs are kept together and the associated Ritz estimates are located in
!     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------
!
!\Data Distribution Note:
!
!  Fortran-D syntax:
!  ================
!  Real resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!  decompose  d1(n), d2(n,ncv)
!  align      resid(i) with d1(i)
!  align      v(i,j)   with d2(i,j)
!  align      workd(i) with d1(i)     range (1:n)
!  align      workd(i) with d1(i-n)   range (n+1:2*n)
!  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
!  distribute d1(block), d2(block,:)
!  replicated workl(lworkl)
!
!  Cray MPP syntax:
!  ===============
!  Real  resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
!  shared     resid(block), v(block,:), workd(block,:)
!  replicated workl(lworkl)
!
!  CM2/CM5 syntax:
!  ==============
!
!-----------------------------------------------------------------------
!
!     include   'ex-nonsym.doc'
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
!     Real Matrices", Linear Algebra and its Applications, vol 88/89,
!     pp 575-595, (1987).
!
!\Routines called:
!     snaup2  ARPACK routine that implements the Implicitly Restarted
!             Arnoldi Iteration.
!     ivout   ARPACK utility routine that prints integers.
!     svout   ARPACK utility routine that prints vectors.
!     slamch  LAPACK routine that determines machine constants.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/16/93: Version '1.1'
!
!\SCCS Information: @(#)
! FILE: naupd.F   SID: 2.10   DATE OF SID: 08/23/02   RELEASE: 2
!
!\Remarks
!
!\EndLib
!
!-----------------------------------------------------------------------
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Real &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(11), ipntr(14)
      Real &
                 resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    bounds, ierr, ih, iq, ishift, iupd, iw, &
                 ldh, ldq, mode, msglvl, mxiter, nb, &
                 nev0, next, np, ritzi, ritzr, j
      save       bounds, ih, iq, ishift, iupd, iw, ldh, ldq, &
                 levec, mode, msglvl, mxiter, nb, nev0, next, &
                 np, ritzi, ritzr
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   snaup2, svout, ivout, sstatn
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slamch
      external   slamch
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido == 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call sstatn
         call cpu_time (t0)
         msglvl = mnaupd
!
!        %----------------%
!        | Error checking |
!        %----------------%
!
         ierr   = 0
         ishift = iparam(1)
!         levec  = iparam(2)
         mxiter = iparam(3)
!         nb     = iparam(4)
         nb     = 1
!
!        %--------------------------------------------%
!        | Revision 2 performs only implicit restart. |
!        %--------------------------------------------%
!
         iupd   = 1
         mode   = iparam(7)
!
         if (n <= 0) then
            ierr = -1
         else if (nev <= 0) then
            ierr = -2
         else if (ncv <= nev+1 .or.  ncv > n) then
            ierr = -3
         else if (mxiter <=          0) then
            ierr = 4
         else if (which /= 'LM' .and. &
             which /= 'SM' .and. &
             which /= 'LR' .and. &
             which /= 'SR' .and. &
             which /= 'LI' .and. &
             which /= 'SI') then
            ierr = -5
         else if (bmat /= 'I' .and. bmat .ne. 'G') then
            ierr = -6
         else if (lworkl < 3*ncv**2 + 6*ncv) then
            ierr = -7
         else if (mode < 1 .or. mode > 4) then
            ierr = -10
         else if (mode == 1 .and. bmat .eq. 'G') then
            ierr = -11
         else if (ishift < 0 .or. ishift > 1) then
            ierr = -12
         end if
!
!        %------------%
!        | Error Exit |
!        %------------%
!
         if (ierr /= 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
!
!        %------------------------%
!        | Set default parameters |
!        %------------------------%
!
         if (nb <= 0)                        nb = 1
         if (tol <= zero)                  tol = slamch('EpsMach')
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        | NEV0 is the local variable designating the   |
!        | size of the invariant subspace desired.      |
!        %----------------------------------------------%
!
         np     = ncv - nev
         nev0   = nev
!
!        %-----------------------------%
!        | Zero out internal workspace |
!        %-----------------------------%
!
         do j = 1, 3*ncv**2 + 6*ncv
            workl(j) = zero
         end do
!
!        %-------------------------------------------------------------%
!        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
!        | etc... and the remaining workspace.                         |
!        | Also update pointer to be used on output.                   |
!        | Memory is laid out as follows:                              |
!        | workl(1:ncv*ncv) := generated Hessenberg matrix             |
!        | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary        |
!        |                                   parts of ritz values      |
!        | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds        |
!        | workl(ncv*ncv+3*ncv+1:2*ncv*ncv+3*ncv) := rotation matrix Q |
!        | workl(2*ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) := workspace       |
!        | The final workspace is needed by subroutine sneigh called   |
!        | by snaup2. Subroutine sneigh calls LAPACK routines for      |
!        | calculating eigenvalues and the last row of the eigenvector |
!        | matrix.                                                     |
!        %-------------------------------------------------------------%
!
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritzr  = ih     + ldh*ncv
         ritzi  = ritzr  + ncv
         bounds = ritzi  + ncv
         iq     = bounds + ncv
         iw     = iq     + ldq*ncv
         next   = iw     + ncv**2 + 3*ncv
!
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritzr
         ipntr(7) = ritzi
         ipntr(8) = bounds
         ipntr(14) = iw
!
      end if
!
!     %-------------------------------------------------------%
!     | Carry out the Implicitly restarted Arnoldi Iteration. |
!     %-------------------------------------------------------%
!
      call snaup2 &
         ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd, &
           ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritzr), &
           workl(ritzi), workl(bounds), workl(iq), ldq, workl(iw), &
           ipntr, workd, info )
!
!     %--------------------------------------------------%
!     | ido /= 99 implies use of reverse communication |
!     | to compute operations involving OP or shifts.    |
!     %--------------------------------------------------%
!
      if (ido == 3) iparam(8) = np
      if (ido /= 99) go to 9000
!
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
!
!     %------------------------------------%
!     | Exit if there was an informational |
!     | error within snaup2.               |
!     %------------------------------------%
!
      if (info < 0) go to 9000
      if (info == 2) info = 3
!
      if (msglvl > 0) then
         call ivout (logfil, 1, mxiter, ndigit, &
                     '_naupd: Number of update iterations taken')
         call ivout (logfil, 1, np, ndigit, &
                     '_naupd: Number of wanted "converged" Ritz values')
         call svout (logfil, np, workl(ritzr), ndigit, &
                     '_naupd: Real part of the final Ritz values')
         call svout (logfil, np, workl(ritzi), ndigit, &
                     '_naupd: Imaginary part of the final Ritz values')
         call svout (logfil, np, workl(bounds), ndigit, &
                     '_naupd: Associated Ritz estimates')
      end if
!
      call cpu_time (t1)
      tnaupd = t1 - t0
!
      if (msglvl > 0) then
!
!        %--------------------------------------------------------%
!        | Version Number & Version Date are defined in version.h |
!        %--------------------------------------------------------%
!
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt, &
                        tmvopx, tmvbx, tnaupd, tnaup2, tnaitr, titref, &
                        tgetv0, tneigh, tngets, tnapps, tnconv, trvec
 1000    format (//, &
            5x, '=============================================',/ &
            5x, '= Nonsymmetric implicit Arnoldi update code =',/ &
            5x, '= Version Number: ', ' 2.4', 21x, ' =',/ &
            5x, '= Version Date:   ', ' 07/31/96', 16x,   ' =',/ &
            5x, '=============================================',/ &
            5x, '= Summary of timing statistics              =',/ &
            5x, '=============================================',//)
 1100    format ( &
            5x, 'Total number update iterations             = ', i5,/ &
            5x, 'Total number of OP*x operations            = ', i5,/ &
            5x, 'Total number of B*x operations             = ', i5,/ &
            5x, 'Total number of reorthogonalization steps  = ', i5,/ &
            5x, 'Total number of iterative refinement steps = ', i5,/ &
            5x, 'Total number of restart steps              = ', i5,/ &
            5x, 'Total time in user OP*x operation          = ', f12.6,/ &
            5x, 'Total time in user B*x operation           = ', f12.6,/ &
            5x, 'Total time in Arnoldi update routine       = ', f12.6,/ &
            5x, 'Total time in naup2 routine                = ', f12.6,/ &
            5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/ &
            5x, 'Total time in reorthogonalization phase    = ', f12.6,/ &
            5x, 'Total time in (re)start vector generation  = ', f12.6,/ &
            5x, 'Total time in Hessenberg eig. subproblem   = ', f12.6,/ &
            5x, 'Total time in getting the shifts           = ', f12.6,/ &
            5x, 'Total time in applying the shifts          = ', f12.6,/ &
            5x, 'Total time in convergence testing          = ', f12.6,/ &
            5x, 'Total time in computing final Ritz vectors = ', f12.6/)
      end if

 9000 continue

      return
      end
      subroutine snconv (n, ritzr, ritzi, bounds, tol, nconv)
!
!! SNCONV does convergence testing for nonsymmetric Arnoldi eigenvalues.
!
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: snconv
!
!\Description:
!  Convergence testing for the nonsymmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call snconv
!     ( N, RITZR, RITZI, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZR,  Real arrays of length N.  (INPUT)
!  RITZI   Real and imaginary parts of the Ritz values to be checked
!          for convergence.

!  BOUNDS  Real array of length N.  (INPUT)
!          Ritz estimates for the Ritz values in RITZR and RITZI.
!
!  TOL     Real scalar.  (INPUT)
!          Desired backward error for a Ritz value to be considered
!          "converged".
!
!  NCONV   Integer scalar.  (OUTPUT)
!          Number of "converged" Ritz values.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     slamch  LAPACK routine that determines machine constants.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: nconv.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     1. xxxx
!
!\EndLib
!
!-----------------------------------------------------------------------
!
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    n, nconv
      Real &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%

      Real &
                 ritzr(n), ritzi(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i
      Real &
                 temp, eps23
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slapy2, slamch
      external   slapy2, slamch

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------------------------------------%
!     | Convergence test: unlike in the symmetric code, I am not    |
!     | using things like refined error bounds and gap condition    |
!     | because I don't know the exact equivalent concept.          |
!     |                                                             |
!     | Instead the i-th Ritz value is considered "converged" when: |
!     |                                                             |
!     |     bounds(i) <= ( TOL * | ritz | )                       |
!     |                                                             |
!     | for some appropriate choice of norm.                        |
!     %-------------------------------------------------------------%
!
      call cpu_time (t0)
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
      eps23 = slamch('Epsilon-Machine')
      eps23 = eps23**(2.0E+0 / 3.0E+0)
!
      nconv  = 0
      do i = 1, n
         temp = max( eps23, slapy2( ritzr(i), ritzi(i) ) )
         if (bounds(i) <= tol*temp)   nconv = nconv + 1
      end do

      call cpu_time (t1)
      tnconv = tnconv + (t1 - t0)

      return
      end
      subroutine sneigh (rnorm, n, h, ldh, ritzr, ritzi, bounds, &
                         q, ldq, workl, ierr)
!
!! SNEIGH computes the eigenvalues of the current upper Hessenberg matrix.
!
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: sneigh
!
!\Description:
!  Compute the eigenvalues of the current upper Hessenberg matrix
!  and the corresponding Ritz estimates given the current residual norm.
!
!\Usage:
!  call sneigh
!     ( RNORM, N, H, LDH, RITZR, RITZI, BOUNDS, Q, LDQ, WORKL, IERR )
!
!\Arguments
!  RNORM   Real scalar.  (INPUT)
!          Residual norm corresponding to the current upper Hessenberg
!          matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the matrix H.
!
!  H       Real N by N array.  (INPUT)
!          H contains the current upper Hessenberg matrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZR,  Real arrays of length N.  (OUTPUT)
!  RITZI   On output, RITZR(1:N) (resp. RITZI(1:N)) contains the real
!          (respectively imaginary) parts of the eigenvalues of H.
!
!  BOUNDS  Real array of length N.  (OUTPUT)
!          On output, BOUNDS contains the Ritz estimates associated with
!          the eigenvalues RITZR and RITZI.  This is equal to RNORM
!          times the last components of the eigenvectors corresponding
!          to the eigenvalues in RITZR and RITZI.
!
!  Q       Real N by N array.  (WORKSPACE)
!          Workspace needed to store the eigenvectors of H.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Real work array of length N**2 + 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  This is needed to keep the full Schur form
!          of H and also in the calculation of the eigenvectors of H.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from slaqrb or strevc.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     slaqrb  ARPACK routine to compute the real Schur form of an
!             upper Hessenberg matrix and last row of the Schur vectors.
!     smout   ARPACK utility routine that prints matrices
!     svout   ARPACK utility routine that prints vectors.
!     slacpy  LAPACK matrix copy routine.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     strevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper quasi-triangular form
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     scopy   Level 1 BLAS that copies one vector to another .
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sscal   Level 1 BLAS that scales a vector.
!
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: neigh.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    ierr, n, ldh, ldq
      Real &
                 rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 bounds(n), h(ldh,n), q(ldq,n), ritzi(n), ritzr(n), &
                 workl(n*(n+3))
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %------------------------%
!     | Local Scalars & Arrays |
!     %------------------------%
!
      logical    select(1)
      integer    i, iconj, msglvl
      Real &
                 temp, vl(1)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy, slacpy, slaqrb, strevc, svout
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slapy2, snrm2
      external   slapy2, snrm2
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic  abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = mneigh
!
      if (msglvl > 2) then
          call smout (logfil, n, n, h, ldh, ndigit, &
               '_neigh: Entering upper Hessenberg matrix H ')
      end if
!
!     %-----------------------------------------------------------%
!     | 1. Compute the eigenvalues, the last components of the    |
!     |    corresponding Schur vectors and the full Schur form T  |
!     |    of the current upper Hessenberg matrix H.              |
!     | slaqrb returns the full Schur form of H in WORKL(1:N**2)  |
!     | and the last components of the Schur vectors in BOUNDS.   |
!     %-----------------------------------------------------------%
!
      call slacpy ('All', n, n, h, ldh, workl, n)
      call slaqrb (.true., n, 1, n, workl, n, ritzr, ritzi, bounds, &
                   ierr)
      if (ierr /= 0) go to 9000
!
      if (msglvl > 1) then
         call svout (logfil, n, bounds, ndigit, &
                    '_neigh: last row of the Schur matrix for H')
      end if
!
!     %-----------------------------------------------------------%
!     | 2. Compute the eigenvectors of the full Schur form T and  |
!     |    apply the last components of the Schur vectors to get  |
!     |    the last components of the corresponding eigenvectors. |
!     | Remember that if the i-th and (i+1)-st eigenvalues are    |
!     | complex conjugate pairs, then the real & imaginary part   |
!     | of the eigenvector components are split across adjacent   |
!     | columns of Q.                                             |
!     %-----------------------------------------------------------%
!
      call strevc ('R', 'A', select, n, workl, n, vl, n, q, ldq, &
                   n, n, workl(n*n+1), ierr)
!
      if (ierr /= 0) go to 9000
!
!     %------------------------------------------------%
!     | Scale the returning eigenvectors so that their |
!     | euclidean norms are all one. LAPACK subroutine |
!     | strevc returns each eigenvector normalized so  |
!     | that the element of largest magnitude has      |
!     | magnitude 1; here the magnitude of a complex   |
!     | number (x,y) is taken to be |x| + |y|.         |
!     %------------------------------------------------%
!
      iconj = 0

      do i=1, n

         if ( abs( ritzi(i) ) <= zero ) then
!
!           %----------------------%
!           | Real eigenvalue case |
!           %----------------------%
!
            temp = snrm2( n, q(1,i), 1 )
            call sscal ( n, one / temp, q(1,i), 1 )
         else
!
!           %-------------------------------------------%
!           | Complex conjugate pair case. Note that    |
!           | since the real and imaginary part of      |
!           | the eigenvector are stored in consecutive |
!           | columns, we further normalize by the      |
!           | square root of two.                       |
!           %-------------------------------------------%
!
            if (iconj == 0) then
               temp = slapy2( snrm2( n, q(1,i), 1 ), &
                              snrm2( n, q(1,i+1), 1 ) )
               call sscal ( n, one / temp, q(1,i), 1 )
               call sscal ( n, one / temp, q(1,i+1), 1 )
               iconj = 1
            else
               iconj = 0
            end if
         end if

      end do

      call sgemv ('T', n, n, one, q, ldq, bounds, 1, zero, workl, 1)

      if (msglvl > 1) then
         call svout (logfil, n, workl, ndigit, &
                    '_neigh: Last row of the eigenvector matrix for H')
      end if
!
!     %----------------------------%
!     | Compute the Ritz estimates |
!     %----------------------------%
!
      iconj = 0

      do i = 1, n

         if ( abs( ritzi(i) ) <= zero ) then
!
!           %----------------------%
!           | Real eigenvalue case |
!           %----------------------%
!
            bounds(i) = rnorm * abs( workl(i) )
         else
!
!           %-------------------------------------------%
!           | Complex conjugate pair case. Note that    |
!           | since the real and imaginary part of      |
!           | the eigenvector are stored in consecutive |
!           | columns, we need to take the magnitude    |
!           | of the last components of the two vectors |
!           %-------------------------------------------%
!
            if (iconj == 0) then
               bounds(i) = rnorm * slapy2( workl(i), workl(i+1) )
               bounds(i+1) = bounds(i)
               iconj = 1
            else
               iconj = 0
            end if
         end if

      end do

      if (msglvl > 2) then
         call svout (logfil, n, ritzr, ndigit, &
                    '_neigh: Real part of the eigenvalues of H')
         call svout (logfil, n, ritzi, ndigit, &
                    '_neigh: Imaginary part of the eigenvalues of H')
         call svout (logfil, n, bounds, ndigit, &
                    '_neigh: Ritz estimates for the eigenvalues of H')
      end if
!
      call cpu_time (t1)
      tneigh = tneigh + (t1 - t0)
!
 9000 continue
      return
!
!     %---------------%
!     | End of sneigh |
!     %---------------%
!
      end
      subroutine sneupd(rvec , howmny, select, dr    , di, &
                         z    , ldz   , sigmar, sigmai, workev, &
                         bmat , n     , which , nev   , tol, &
                         resid, ncv   , v     , ldv   , iparam, &
                         ipntr, workd , workl , lworkl, info)
!
!! SNEUPD returns the converged approximate eigenvalues.
!
!\BeginDoc
!
!\Name: sneupd
!
!\Description:
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) The corresponding approximate eigenvectors;
!
!      (2) An orthonormal basis for the associated approximate
!          invariant subspace;
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  basis is always computed.  There is an additional storage cost of n*nev
!  if both are requested (in this case a separate array Z must be supplied).
!
!  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
!  are derived from approximate eigenvalues and eigenvectors of
!  of the linear operator OP prescribed by the MODE selection in the
!  call to SNAUPD.  SNAUPD must be called before this routine is called.
!  These approximate eigenvalues and vectors are commonly called Ritz
!  values and Ritz vectors respectively.  They are referred to as such
!  in the comments that follow.  The computed orthonormal basis for the
!  invariant subspace corresponding to these Ritz values is referred to as a
!  Schur basis.
!
!  See documentation in the header of the subroutine SNAUPD for
!  definition of OP as well as other terms and the relation of computed
!  Ritz values and Ritz vectors of OP with respect to the given problem
!  A*z = lambda*B*z.  For a brief description, see definitions of
!  IPARAM(7), MODE and WHICH in the documentation of SNAUPD.
!
!\Usage:
!  call sneupd
!     ( RVEC, HOWMNY, SELECT, DR, DI, Z, LDZ, SIGMAR, SIGMAI, WORKEV, BMAT,
!       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL,
!       LWORKL, INFO )
!
!\Arguments:
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether a basis for the invariant subspace corresponding
!          to the converged Ritz value approximations for the eigenproblem
!          A*z = lambda*B*z is computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute the Ritz vectors or Schur vectors.
!                                See Remarks below.
!
!  HOWMNY  Character*1  (INPUT)
!          Specifies the form of the basis for the invariant subspace
!          corresponding to the converged Ritz values that is to be computed.
!
!          = 'A': Compute NEV Ritz vectors;
!          = 'P': Compute NEV Schur vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value (DR(j), DI(j)), SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A' or 'P', SELECT is used as internal workspace.
!
!  DR      Real  array of dimension NEV+1.  (OUTPUT)
!          If IPARAM(7) = 1,2 or 3 and SIGMAI=0.0  then on exit: DR contains
!          the real part of the Ritz  approximations to the eigenvalues of
!          A*z = lambda*B*z.
!          If IPARAM(7) = 3, 4 and SIGMAI is not equal to zero, then on exit:
!          DR contains the real part of the Ritz values of OP computed by
!          SNAUPD. A further computation must be performed by the user
!          to transform the Ritz values computed for OP by SNAUPD to those
!          of the original system A*z = lambda*B*z. See remark 3 below.
!
!  DI      Real  array of dimension NEV+1.  (OUTPUT)
!          On exit, DI contains the imaginary part of the Ritz value
!          approximations to the eigenvalues of A*z = lambda*B*z associated
!          with DR.
!
!          NOTE: When Ritz values are complex, they will come in complex
!                conjugate pairs.  If eigenvectors are requested, the
!                corresponding Ritz vectors will also come in conjugate
!                pairs and the real and imaginary parts of these are
!                represented in two consecutive columns of the array Z
!                (see below).
!
!  Z       Real  N by NEV+1 array if RVEC = .TRUE. and HOWMNY = 'A'. (OUTPUT)
!          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
!          Z represent approximate eigenvectors (Ritz vectors) corresponding
!          to the NCONV=IPARAM(5) Ritz values for eigensystem
!          A*z = lambda*B*z.
!
!          The complex Ritz vector associated with the Ritz value
!          with positive imaginary part is stored in two consecutive
!          columns.  The first column holds the real part of the Ritz
!          vector and the second column holds the imaginary part.  The
!          Ritz vector associated with the Ritz value with negative
!          imaginary part is simply the complex conjugate of the Ritz vector
!          associated with the positive imaginary part.
!
!          If  RVEC = .FALSE. or HOWMNY = 'P', then Z is not referenced.
!
!          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
!          the array Z may be set equal to first NEV+1 columns of the Arnoldi
!          basis array V computed by SNAUPD.  In this case the Arnoldi basis
!          will be destroyed and overwritten with the eigenvector basis.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ >= max( 1, N ).  In any case,  LDZ >= 1.
!
!  SIGMAR  Real   (INPUT)
!          If IPARAM(7) = 3 or 4, represents the real part of the shift.
!          Not referenced if IPARAM(7) = 1 or 2.
!
!  SIGMAI  Real   (INPUT)
!          If IPARAM(7) = 3 or 4, represents the imaginary part of the shift.
!          Not referenced if IPARAM(7) = 1 or 2. See remark 3 below.
!
!  WORKEV  Real  work array of dimension 3*NCV.  (WORKSPACE)
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to SNAUPD that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to SNEUPD following the last call
!         to SNAUPD.  These arguments MUST NOT BE MODIFIED between
!         the the last call to SNAUPD and the call to SNEUPD.
!
!  Three of these parameters (V, WORKL, INFO) are also output parameters:
!
!  V       Real  N by NCV array.  (INPUT/OUTPUT)
!
!          Upon INPUT: the NCV columns of V contain the Arnoldi basis
!                      vectors for OP as constructed by SNAUPD .
!
!          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
!                       contain approximate Schur vectors that span the
!                       desired invariant subspace.  See Remark 2 below.
!
!          NOTE: If the array Z has been set equal to first NEV+1 columns
!          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
!          Arnoldi basis held by V has been overwritten by the desired
!          Ritz vectors.  If a separate array Z has been passed then
!          the first NCONV=IPARAM(5) columns of V will contain approximate
!          Schur vectors that span the desired invariant subspace.
!
!  WORKL   Real  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:ncv*ncv+3*ncv) contains information obtained in
!          snaupd.  They are not changed by sneupd.
!          WORKL(ncv*ncv+3*ncv+1:3*ncv*ncv+6*ncv) holds the
!          real and imaginary part of the untransformed Ritz values,
!          the upper quasi-triangular matrix for H, and the
!          associated matrix representation of the invariant subspace for H.
!
!          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
!          of the above information computed by sneupd.
!          -------------------------------------------------------------
!          IPNTR(9):  pointer to the real part of the NCV RITZ values of the
!                     original system.
!          IPNTR(10): pointer to the imaginary part of the NCV RITZ values of
!                     the original system.
!          IPNTR(11): pointer to the NCV corresponding error bounds.
!          IPNTR(12): pointer to the NCV by NCV upper quasi-triangular
!                     Schur matrix for H.
!          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
!                     of the upper Hessenberg matrix H. Only referenced by
!                     sneupd if RVEC = .TRUE. See Remark 2 below.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!
!          =  0: Normal exit.
!
!          =  1: The Schur form computed by LAPACK routine slahqr
!                could not be reordered by LAPACK routine strsen.
!                Re-enter subroutine sneupd with IPARAM(5)=NCV and
!                increase the size of the arrays DR and DI to have
!                dimension at least dimension NCV and allocate at least NCV
!                columns for Z. NOTE: Not necessary if Z and V share
!                the same space. Please notify the authors if this error
!                occurs.
!
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV-NEV >= 2 and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from calculation of a real Schur form.
!                Informational error from LAPACK routine slahqr.
!          = -9: Error return from calculation of eigenvectors.
!                Informational error from LAPACK routine strevc.
!          = -10: IPARAM(7) must be 1,2,3,4.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: HOWMNY = 'S' not yet implemented
!          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
!          = -14: SNAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: DNEUPD got a different count of the number of converged
!                 Ritz values than DNAUPD got.  This indicates the user
!                 probably made an error in passing data from DNAUPD to
!                 DNEUPD or that the data was modified before entering
!                 DNEUPD
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett & Y. Saad, "Complex Shift and Invert Strategies for
!     Real Matrices", Linear Algebra and its Applications, vol 88/89,
!     pp 575-595, (1987).
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers.
!     smout   ARPACK utility routine that prints matrices
!     svout   ARPACK utility routine that prints vectors.
!     sgeqr2  LAPACK routine that computes the QR factorization of
!             a matrix.
!     slacpy  LAPACK matrix copy routine.
!     slahqr  LAPACK routine to compute the real Schur form of an
!             upper Hessenberg matrix.
!     slamch  LAPACK routine that determines machine constants.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     slaset  LAPACK matrix initialization routine.
!     sorm2r  LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     strevc  LAPACK routine to compute the eigenvectors of a matrix
!             in upper quasi-triangular form.
!     strsen  LAPACK routine that re-orders the Schur form.
!     strmm   Level 3 BLAS matrix times an upper triangular matrix.
!     sger    Level 2 BLAS rank one update to a matrix.
!     scopy   Level 1 BLAS that copies one vector to another .
!     sdot    Level 1 BLAS that computes the scalar product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sscal   Level 1 BLAS that scales a vector.
!
!\Remarks
!
!  1. Currently only HOWMNY = 'A' and 'P' are implemented.
!
!     Let trans(X) denote the transpose of X.
!
!  2. Schur vectors are an orthogonal representation for the basis of
!     Ritz vectors. Thus, their numerical properties are often superior.
!     If RVEC = .TRUE. then the relationship
!             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
!     trans(V(:,1:IPARAM(5))) * V(:,1:IPARAM(5)) = I are approximately
!     satisfied. Here T is the leading submatrix of order IPARAM(5) of the
!     real upper quasi-triangular matrix stored workl(ipntr(12)). That is,
!     T is block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!     each 2-by-2 diagonal block has its diagonal elements equal and its
!     off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
!     diagonal block is a complex conjugate pair of Ritz values. The real
!     Ritz values are stored on the diagonal of T.
!
!  3. If IPARAM(7) = 3 or 4 and SIGMAI is not equal zero, then the user must
!     form the IPARAM(5) Rayleigh quotients in order to transform the Ritz
!     values computed by SNAUPD for OP to those of A*z = lambda*B*z.
!     Set RVEC = .true. and HOWMNY = 'A', and
!     compute
!           trans(Z(:,I)) * A * Z(:,I) if DI(I) = 0.
!     If DI(I) is not equal to zero and DI(I+1) = - D(I),
!     then the desired real and imaginary parts of the Ritz value are
!           trans(Z(:,I)) * A * Z(:,I) +  trans(Z(:,I+1)) * A * Z(:,I+1),
!           trans(Z(:,I)) * A * Z(:,I+1) -  trans(Z(:,I+1)) * A * Z(:,I),
!     respectively.
!     Another possibility is to set RVEC = .true. and HOWMNY = 'P' and
!     compute trans(V(:,1:IPARAM(5))) * A * V(:,1:IPARAM(5)) and then an upper
!     quasi-triangular matrix of order IPARAM(5) is computed. See remark
!     2 above.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: neupd.F   SID: 2.7   DATE OF SID: 09/20/00   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Real &
                 sigmar, sigmai, tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(11), ipntr(14)
      logical    select(ncv)
      Real &
                 dr(nev+1)    , di(nev+1), resid(n)  , &
                 v(ldv,ncv)   , z(ldz,*) , workd(3*n), &
                 workl(lworkl), workev(3*ncv)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0 , zero = 0.0E+0 )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  type*6
      integer    bounds, ierr  , ih    , ihbds   , &
                 iheigr, iheigi, iconj , nconv   , &
                 invsub, iuptri, iwev  , iwork(1), &
                 j     , k     , ldh   , ldq     , &
                 mode  , msglvl, outncv, ritzr   , &
                 ritzi , wri   , wrr   , irr     , &
                 iri   , ibd   , ishift, numcnv  , &
                 np    , jj
      logical    reord
      Real &
                 conds  , rnorm, sep  , temp, &
                 vl(1,1), temp1, eps23
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy , sger  , sgeqr2, slacpy, &
                 slahqr, slaset, smout , sorm2r, &
                 strevc, strmm , strsen, sscal , &
                 svout , ivout
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slapy2, snrm2, slamch, sdot
      external   slapy2, snrm2, slamch, sdot
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs, min, sqrt
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
      msglvl = mneupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
!
!     %---------------------------------%
!     | Get machine dependent constant. |
!     %---------------------------------%
!
      eps23 = slamch('Epsilon-Machine')
      eps23 = eps23**(2.0E+0  / 3.0E+0 )
!
!     %--------------%
!     | Quick return |
!     %--------------%
!
      ierr = 0
!
      if (nconv <= 0) then
         ierr = -14
      else if (n <= 0) then
         ierr = -1
      else if (nev <= 0) then
         ierr = -2
      else if (ncv <= nev+1 .or.  ncv > n) then
         ierr = -3
      else if (which /= 'LM' .and. &
              which /= 'SM' .and. &
              which /= 'LR' .and. &
              which /= 'SR' .and. &
              which /= 'LI' .and. &
              which /= 'SI') then
         ierr = -5
      else if (bmat /= 'I' .and. bmat .ne. 'G') then
         ierr = -6
      else if (lworkl < 3*ncv**2 + 6*ncv) then
         ierr = -7
      else if ( (howmny /= 'A' .and. &
                 howmny /= 'P' .and. &
                 howmny /= 'S') .and. rvec ) then
         ierr = -13
      else if (howmny == 'S' ) then
         ierr = -12
      end if
!
      if (mode == 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode == 3 .and. sigmai .eq. zero) then
         type = 'SHIFTI'
      else if (mode == 3 ) then
         type = 'REALPT'
      else if (mode == 4 ) then
         type = 'IMAGPT'
      else
                                              ierr = -10
      end if
      if (mode == 1 .and. bmat .eq. 'G')    ierr = -11
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      if (ierr /= 0) then
         info = ierr
         go to 9000
      end if
!
!     %--------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q   |
!     | etc... and the remaining workspace.                    |
!     | Also update pointer to be used on output.              |
!     | Memory is laid out as follows:                         |
!     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
!     | workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary   |
!     |                                   parts of ritz values |
!     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds   |
!     %--------------------------------------------------------%
!
!     %-----------------------------------------------------------%
!     | The following is used and set by SNEUPD.                  |
!     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
!     |                             real part of the Ritz values. |
!     | workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed |
!     |                        imaginary part of the Ritz values. |
!     | workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed |
!     |                           error bounds of the Ritz values |
!     | workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper |
!     |                             quasi-triangular matrix for H |
!     | workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the    |
!     |       associated matrix representation of the invariant   |
!     |       subspace for H.                                     |
!     | GRAND total of NCV * ( 3 * NCV + 6 ) locations.           |
!     %-----------------------------------------------------------%
!
      ih     = ipntr(5)
      ritzr  = ipntr(6)
      ritzi  = ipntr(7)
      bounds = ipntr(8)
      ldh    = ncv
      ldq    = ncv
      iheigr = bounds + ldh
      iheigi = iheigr + ldh
      ihbds  = iheigi + ldh
      iuptri = ihbds  + ldh
      invsub = iuptri + ldh*ncv
      ipntr(9)  = iheigr
      ipntr(10) = iheigi
      ipntr(11) = ihbds
      ipntr(12) = iuptri
      ipntr(13) = invsub
      wrr = 1
      wri = ncv + 1
      iwev = wri + ncv
!
!     %-----------------------------------------%
!     | irr points to the REAL part of the Ritz |
!     |     values computed by _neigh before    |
!     |     exiting _naup2.                     |
!     | iri points to the IMAGINARY part of the |
!     |     Ritz values computed by _neigh      |
!     |     before exiting _naup2.              |
!     | ibd points to the Ritz estimates        |
!     |     computed by _neigh before exiting   |
!     |     _naup2.                             |
!     %-----------------------------------------%
!
      irr = ipntr(14)+ncv*ncv
      iri = irr+ncv
      ibd = iri+ncv
!
!     %------------------------------------%
!     | RNORM is B-norm of the RESID(1:N). |
!     %------------------------------------%
!
      rnorm = workl(ih+2)
      workl(ih+2) = zero

      if (msglvl > 2) then
         call svout(logfil, ncv, workl(irr), ndigit, &
         '_neupd: Real part of Ritz values passed in from _NAUPD.')
         call svout(logfil, ncv, workl(iri), ndigit, &
         '_neupd: Imag part of Ritz values passed in from _NAUPD.')
         call svout(logfil, ncv, workl(ibd), ndigit, &
         '_neupd: Ritz estimates passed in from _NAUPD.')
      end if

      if (rvec) then

         reord = .false.
!
!        %---------------------------------------------------%
!        | Use the temporary bounds array to store indices   |
!        | These will be used to mark the select array later |
!        %---------------------------------------------------%
!
         do j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
         end do
!
!        %-------------------------------------%
!        | Select the wanted Ritz values.      |
!        | Sort the Ritz values so that the    |
!        | wanted ones appear at the tailing   |
!        | NEV positions of workl(irr) and     |
!        | workl(iri).  Move the corresponding |
!        | error estimates in workl(bound)     |
!        | accordingly.                        |
!        %-------------------------------------%
!
         np     = ncv - nev
         ishift = 0
         call sngets(ishift       , which     , nev       , &
                      np           , workl(irr), workl(iri), &
                      workl(bounds), workl     , workl(np+1))
!
         if (msglvl > 2) then
            call svout(logfil, ncv, workl(irr), ndigit, &
            '_neupd: Real part of Ritz values after calling _NGETS.')
            call svout(logfil, ncv, workl(iri), ndigit, &
            '_neupd: Imag part of Ritz values after calling _NGETS.')
            call svout(logfil, ncv, workl(bounds), ndigit, &
            '_neupd: Ritz value indices after calling _NGETS.')
         end if
!
!        %-----------------------------------------------------%
!        | Record indices of the converged wanted Ritz values  |
!        | Mark the select array for possible reordering       |
!        %-----------------------------------------------------%
!
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23, &
                       slapy2( workl(irr+ncv-j), workl(iri+ncv-j) ))
            jj = workl(bounds + ncv - j)
            if (numcnv < nconv .and. &
                workl(ibd+jj-1) <= tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj > nev) reord = .true.
            endif
   11    continue
!
!        %-----------------------------------------------------------%
!        | Check the count (numcnv) of converged Ritz values with    |
!        | the number (nconv) reported by dnaupd.  If these two      |
!        | are different then there has probably been an error       |
!        | caused by incorrect passing of the dnaupd data.           |
!        %-----------------------------------------------------------%
!
         if (msglvl > 2) then
             call ivout(logfil, 1, numcnv, ndigit, &
                  '_neupd: Number of specified eigenvalues')
             call ivout(logfil, 1, nconv, ndigit, &
                  '_neupd: Number of "converged" eigenvalues')
         end if
!
         if (numcnv /= nconv) then
            info = -15
            go to 9000
         end if
!
!        %-----------------------------------------------------------%
!        | Call LAPACK routine slahqr to compute the real Schur form |
!        | of the upper Hessenberg matrix returned by SNAUPD.        |
!        | Make a copy of the upper Hessenberg matrix.               |
!        | Initialize the Schur vector matrix Q to the identity.     |
!        %-----------------------------------------------------------%
!
         call scopy(ldh*ncv, workl(ih), 1, workl(iuptri), 1)
         call slaset('All', ncv, ncv, &
                      zero , one, workl(invsub), &
                      ldq)
         call slahqr(.true., .true.       , ncv, &
                      1     , ncv          , workl(iuptri), &
                      ldh   , workl(iheigr), workl(iheigi), &
                      1     , ncv          , workl(invsub), &
                      ldq   , ierr)
         call scopy(ncv         , workl(invsub+ncv-1), ldq, &
                     workl(ihbds), 1)
!
         if (ierr /= 0) then
            info = -8
            go to 9000
         end if
!
         if (msglvl > 1) then
            call svout(logfil, ncv, workl(iheigr), ndigit, &
                 '_neupd: Real part of the eigenvalues of H')
            call svout(logfil, ncv, workl(iheigi), ndigit, &
                 '_neupd: Imaginary part of the Eigenvalues of H')
            call svout(logfil, ncv, workl(ihbds), ndigit, &
                 '_neupd: Last row of the Schur vector matrix')
            if (msglvl > 3) then
               call smout(logfil       , ncv, ncv   , &
                           workl(iuptri), ldh, ndigit, &
                    '_neupd: The upper quasi-triangular matrix ')
            end if
         end if
!
         if (reord) then
!
!           %-----------------------------------------------------%
!           | Reorder the computed upper quasi-triangular matrix. |
!           %-----------------------------------------------------%
!
            call strsen('None'       , 'V'          , &
                         select       , ncv          , &
                         workl(iuptri), ldh          , &
                         workl(invsub), ldq          , &
                         workl(iheigr), workl(iheigi), &
                         nconv        , conds        , &
                         sep          , workl(ihbds) , &
                         ncv          , iwork        , &
                         1            , ierr)

            if (ierr == 1) then
               info = 1
               go to 9000
            end if
!
            if (msglvl > 2) then
                call svout(logfil, ncv, workl(iheigr), ndigit, &
                 '_neupd: Real part of the eigenvalues of H--reordered')
                call svout(logfil, ncv, workl(iheigi), ndigit, &
                 '_neupd: Imag part of the eigenvalues of H--reordered')
                if (msglvl > 3) then
                   call smout(logfil       , ncv, ncv   , &
                               workl(iuptri), ldq, ndigit, &
                   '_neupd: Quasi-triangular matrix after re-ordering')
                end if
            end if
!
         end if
!
!        %---------------------------------------%
!        | Copy the last row of the Schur vector |
!        | into workl(ihbds).  This will be used |
!        | to compute the Ritz estimates of      |
!        | converged Ritz values.                |
!        %---------------------------------------%
!
         call scopy(ncv, workl(invsub+ncv-1), ldq, workl(ihbds), 1)
!
!        %----------------------------------------------------%
!        | Place the computed eigenvalues of H into DR and DI |
!        | if a spectral transformation was not used.         |
!        %----------------------------------------------------%
!
         if (type == 'REGULR') then
            call scopy(nconv, workl(iheigr), 1, dr, 1)
            call scopy(nconv, workl(iheigi), 1, di, 1)
         end if
!
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(invsub,ldq).                            |
!        %----------------------------------------------------------%
!
         call sgeqr2(ncv, nconv , workl(invsub), &
                     ldq, workev, workev(ncv+1), &
                     ierr)
!
!        %---------------------------------------------------------%
!        | * Postmultiply V by Q using sorm2r.                     |
!        | * Copy the first NCONV columns of VQ into Z.            |
!        | * Postmultiply Z by R.                                  |
!        | The N by NCONV matrix Z is now a matrix representation  |
!        | of the approximate invariant subspace associated with   |
!        | the Ritz values in workl(iheigr) and workl(iheigi)      |
!        | The first NCONV columns of V are now approximate Schur  |
!        | vectors associated with the real upper quasi-triangular |
!        | matrix of order NCONV in workl(iuptri)                  |
!        %---------------------------------------------------------%
!
         call sorm2r('Right', 'Notranspose', n            , &
                      ncv   , nconv        , workl(invsub), &
                      ldq   , workev       , v            , &
                      ldv   , workd(n+1)   , ierr)
         call slacpy('All', n, nconv, v, ldv, z, ldz)
!
         do j=1, nconv
!
!           %---------------------------------------------------%
!           | Perform both a column and row scaling if the      |
!           | diagonal element of workl(invsub,ldq) is negative |
!           | I'm lazy and don't take advantage of the upper    |
!           | quasi-triangular form of workl(iuptri,ldq)        |
!           | Note that since Q is orthogonal, R is a diagonal  |
!           | matrix consisting of plus or minus ones           |
!           %---------------------------------------------------%
!
            if (workl(invsub+(j-1)*ldq+j-1) < zero) then
               call sscal(nconv, -one, workl(iuptri+j-1), ldq)
               call sscal(nconv, -one, workl(iuptri+(j-1)*ldq), 1)
            end if

         end do
!
         if (howmny == 'A') then
!
!           %--------------------------------------------%
!           | Compute the NCONV wanted eigenvectors of T |
!           | located in workl(iuptri,ldq).              |
!           %--------------------------------------------%
!
            do 30 j=1, ncv
               if (j <= nconv) then
                  select(j) = .true.
               else
                  select(j) = .false.
               end if
 30         continue
!
            call strevc('Right', 'Select'     , select       , &
                         ncv    , workl(iuptri), ldq          , &
                         vl     , 1            , workl(invsub), &
                         ldq    , ncv          , outncv       , &
                         workev , ierr)
!
            if (ierr /= 0) then
                info = -9
                go to 9000
            end if
!
!           %------------------------------------------------%
!           | Scale the returning eigenvectors so that their |
!           | Euclidean norms are all one. LAPACK subroutine |
!           | strevc returns each eigenvector normalized so  |
!           | that the element of largest magnitude has      |
!           | magnitude 1;                                   |
!           %------------------------------------------------%
!
            iconj = 0
            do 40 j=1, nconv
!
               if ( workl(iheigi+j-1) == zero ) then
!
!                 %----------------------%
!                 | real eigenvalue case |
!                 %----------------------%
!
                  temp = snrm2( ncv, workl(invsub+(j-1)*ldq), 1 )
                  call sscal( ncv, one / temp, &
                       workl(invsub+(j-1)*ldq), 1 )
!
               else
!
!                 %-------------------------------------------%
!                 | Complex conjugate pair case. Note that    |
!                 | since the real and imaginary part of      |
!                 | the eigenvector are stored in consecutive |
!                 | columns, we further normalize by the      |
!                 | square root of two.                       |
!                 %-------------------------------------------%
!
                  if (iconj == 0) then
                     temp = slapy2(snrm2(ncv, &
                                         workl(invsub+(j-1)*ldq), &
                                         1), &
                                   snrm2(ncv, &
                                         workl(invsub+j*ldq), &
                                         1))
                     call sscal(ncv, one/temp, &
                                 workl(invsub+(j-1)*ldq), 1 )
                     call sscal(ncv, one/temp, &
                                 workl(invsub+j*ldq), 1 )
                     iconj = 1
                  else
                     iconj = 0
                  end if
!
               end if
!
 40         continue
!
            call sgemv('T', ncv, nconv, one, workl(invsub), &
                       ldq, workl(ihbds), 1, zero,  workev, 1)
!
            iconj = 0
            do 45 j=1, nconv
               if (workl(iheigi+j-1) /= zero) then
!
!                 %-------------------------------------------%
!                 | Complex conjugate pair case. Note that    |
!                 | since the real and imaginary part of      |
!                 | the eigenvector are stored in consecutive |
!                 %-------------------------------------------%
!
                  if (iconj == 0) then
                     workev(j) = slapy2(workev(j), workev(j+1))
                     workev(j+1) = workev(j)
                     iconj = 1
                  else
                     iconj = 0
                  end if
               end if
 45         continue
!
            if (msglvl > 2) then
               call scopy(ncv, workl(invsub+ncv-1), ldq, &
                          workl(ihbds), 1)
               call svout(logfil, ncv, workl(ihbds), ndigit, &
                    '_neupd: Last row of the eigenvector matrix for T')
               if (msglvl > 3) then
                  call smout(logfil, ncv, ncv, workl(invsub), ldq, &
                       ndigit, '_neupd: The eigenvector matrix for T')
               end if
            end if
!
!           %---------------------------------------%
!           | Copy Ritz estimates into workl(ihbds) |
!           %---------------------------------------%
!
            call scopy(nconv, workev, 1, workl(ihbds), 1)
!
!           %---------------------------------------------------------%
!           | Compute the QR factorization of the eigenvector matrix  |
!           | associated with leading portion of T in the first NCONV |
!           | columns of workl(invsub,ldq).                           |
!           %---------------------------------------------------------%
!
            call sgeqr2(ncv, nconv , workl(invsub), &
                         ldq, workev, workev(ncv+1), &
                         ierr)
!
!           %----------------------------------------------%
!           | * Postmultiply Z by Q.                       |
!           | * Postmultiply Z by R.                       |
!           | The N by NCONV matrix Z is now contains the  |
!           | Ritz vectors associated with the Ritz values |
!           | in workl(iheigr) and workl(iheigi).          |
!           %----------------------------------------------%
!
            call sorm2r('Right', 'Notranspose', n            , &
                         ncv  , nconv        , workl(invsub), &
                         ldq  , workev       , z            , &
                         ldz  , workd(n+1)   , ierr)
!
            call strmm('Right'   , 'Upper'       , 'No transpose', &
                        'Non-unit', n            , nconv         , &
                        one       , workl(invsub), ldq           , &
                        z         , ldz)
!
         end if
!
      else
!
!        %------------------------------------------------------%
!        | An approximate invariant subspace is not needed.     |
!        | Place the Ritz values computed SNAUPD into DR and DI |
!        %------------------------------------------------------%
!
         call scopy(nconv, workl(ritzr), 1, dr, 1)
         call scopy(nconv, workl(ritzi), 1, di, 1)
         call scopy(nconv, workl(ritzr), 1, workl(iheigr), 1)
         call scopy(nconv, workl(ritzi), 1, workl(iheigi), 1)
         call scopy(nconv, workl(bounds), 1, workl(ihbds), 1)
      end if
!
!     %------------------------------------------------%
!     | Transform the Ritz values and possibly vectors |
!     | and corresponding error bounds of OP to those  |
!     | of A*x = lambda*B*x.                           |
!     %------------------------------------------------%
!
      if (type == 'REGULR') then
!
         if (rvec) &
            call sscal(ncv, rnorm, workl(ihbds), 1)
!
      else
!
!        %---------------------------------------%
!        |   A spectral transformation was used. |
!        | * Determine the Ritz estimates of the |
!        |   Ritz values in the original system. |
!        %---------------------------------------%
!
         if (type == 'SHIFTI') then
!
            if (rvec) &
               call sscal(ncv, rnorm, workl(ihbds), 1)
!
            do 50 k=1, ncv
               temp = slapy2( workl(iheigr+k-1), &
                              workl(iheigi+k-1) )
               workl(ihbds+k-1) = abs( workl(ihbds+k-1) ) &
                                / temp / temp
 50         continue
!
         else if (type == 'REALPT') then
!
            do 60 k=1, ncv
 60         continue
!
         else if (type == 'IMAGPT') then
!
            do 70 k=1, ncv
 70         continue
!
         end if
!
!        %-----------------------------------------------------------%
!        | *  Transform the Ritz values back to the original system. |
!        |    For TYPE = 'SHIFTI' the transformation is              |
!        |             lambda = 1/theta + sigma                      |
!        |    For TYPE = 'REALPT' or 'IMAGPT' the user must from     |
!        |    Rayleigh quotients or a projection. See remark 3 above.|
!        | NOTES:                                                    |
!        | *The Ritz vectors are not affected by the transformation. |
!        %-----------------------------------------------------------%
!
         if (type == 'SHIFTI') then
!
            do 80 k=1, ncv
               temp = slapy2( workl(iheigr+k-1), &
                              workl(iheigi+k-1) )
               workl(iheigr+k-1) = workl(iheigr+k-1)/temp/temp &
                                 + sigmar
               workl(iheigi+k-1) = -workl(iheigi+k-1)/temp/temp &
                                 + sigmai
 80         continue
!
            call scopy(nconv, workl(iheigr), 1, dr, 1)
            call scopy(nconv, workl(iheigi), 1, di, 1)
!
         else if (type == 'REALPT' .or. type .eq. 'IMAGPT') then
!
            call scopy(nconv, workl(iheigr), 1, dr, 1)
            call scopy(nconv, workl(iheigi), 1, di, 1)
!
         end if
!
      end if
!
      if (type == 'SHIFTI' .and. msglvl > 1) then
         call svout(logfil, nconv, dr, ndigit, &
         '_neupd: Untransformed real part of the Ritz valuess.')
         call svout (logfil, nconv, di, ndigit, &
         '_neupd: Untransformed imag part of the Ritz valuess.')
         call svout(logfil, nconv, workl(ihbds), ndigit, &
         '_neupd: Ritz estimates of untransformed Ritz values.')
      else if (type == 'REGULR' .and. msglvl > 1) then
         call svout(logfil, nconv, dr, ndigit, &
         '_neupd: Real parts of converged Ritz values.')
         call svout (logfil, nconv, di, ndigit, &
         '_neupd: Imag parts of converged Ritz values.')
         call svout(logfil, nconv, workl(ihbds), ndigit, &
         '_neupd: Associated Ritz estimates.')
      end if
!
!     %-------------------------------------------------%
!     | Eigenvector Purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 2.                                   |
!     %-------------------------------------------------%
!
      if (rvec .and. howmny == 'A' .and. type .eq. 'SHIFTI') then
!
!        %------------------------------------------------%
!        | Purify the computed Ritz vectors by adding a   |
!        | little bit of the residual vector:             |
!        |                      T                         |
!        |          resid(:)*( e    s ) / theta           |
!        |                      NCV                       |
!        | where H s = s theta. Remember that when theta  |
!        | has nonzero imaginary part, the corresponding  |
!        | Ritz vector is stored across two columns of Z. |
!        %------------------------------------------------%
!
         iconj = 0
         do 110 j=1, nconv
            if (workl(iheigi+j-1) == zero) then
               workev(j) =  workl(invsub+(j-1)*ldq+ncv-1) / &
                            workl(iheigr+j-1)
            else if (iconj == 0) then
               temp = slapy2( workl(iheigr+j-1), workl(iheigi+j-1) )
               workev(j) = ( workl(invsub+(j-1)*ldq+ncv-1) * &
                             workl(iheigr+j-1) + &
                             workl(invsub+j*ldq+ncv-1) * &
                             workl(iheigi+j-1) ) / temp / temp
               workev(j+1) = ( workl(invsub+j*ldq+ncv-1) * &
                               workl(iheigr+j-1) - &
                               workl(invsub+(j-1)*ldq+ncv-1) * &
                               workl(iheigi+j-1) ) / temp / temp
               iconj = 1
            else
               iconj = 0
            end if
 110     continue
!
!        %---------------------------------------%
!        | Perform a rank one update to Z and    |
!        | purify all the Ritz vectors together. |
!        %---------------------------------------%
!
         call sger(n, nconv, one, resid, 1, workev, 1, z, ldz)
!
      end if
!
 9000 continue
!
      return
!
!     %---------------%
!     | End of SNEUPD |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: sngets
!
!\Description:
!  Given the eigenvalues of the upper Hessenberg matrix H,
!  computes the NP shifts AMU that are zeros of the polynomial of
!  degree NP which filters out components of the unwanted eigenvectors
!  corresponding to the AMU's based on some given criteria.
!
!  NOTE: call this even in the case of user specified shifts in order
!  to sort the eigenvalues, and error bounds of H for later use.
!
!\Usage:
!  call sngets
!     ( ISHIFT, WHICH, KEV, NP, RITZR, RITZI, BOUNDS, SHIFTR, SHIFTI )
!
!\Arguments
!  ISHIFT  Integer.  (INPUT)
!          Method for selecting the implicit shifts at each iteration.
!          ISHIFT = 0: user specified shifts
!          ISHIFT = 1: exact shift with respect to the matrix H.
!
!  WHICH   Character*2.  (INPUT)
!          Shift selection criteria.
!          'LM' -> want the KEV eigenvalues of largest magnitude.
!          'SM' -> want the KEV eigenvalues of smallest magnitude.
!          'LR' -> want the KEV eigenvalues of largest real part.
!          'SR' -> want the KEV eigenvalues of smallest real part.
!          'LI' -> want the KEV eigenvalues of largest imaginary part.
!          'SI' -> want the KEV eigenvalues of smallest imaginary part.
!
!  KEV      Integer.  (INPUT/OUTPUT)
!           INPUT: KEV+NP is the size of the matrix H.
!           OUTPUT: Possibly increases KEV by one to keep complex conjugate
!           pairs together.
!
!  NP       Integer.  (INPUT/OUTPUT)
!           Number of implicit shifts to be computed.
!           OUTPUT: Possibly decreases NP by one to keep complex conjugate
!           pairs together.
!
!  RITZR,  Real array of length KEV+NP.  (INPUT/OUTPUT)
!  RITZI   On INPUT, RITZR and RITZI contain the real and imaginary
!          parts of the eigenvalues of H.
!          On OUTPUT, RITZR and RITZI are sorted so that the unwanted
!          eigenvalues are in the first NP locations and the wanted
!          portion is in the last KEV locations.  When exact shifts are
!          selected, the unwanted part corresponds to the shifts to
!          be applied. Also, if ISHIFT == 1, the unwanted eigenvalues
!          are further sorted so that the ones with largest Ritz values
!          are first.
!
!  BOUNDS  Real array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!  SHIFTR, SHIFTI  *** USE deprecated as of version 2.1. ***
!
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     ssortc  ARPACK sorting routine.
!     scopy   Level 1 BLAS that copies one vector to another .
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: ngets.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\Remarks
!     1. xxxx
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine sngets ( ishift, which, kev, np, ritzr, ritzi, bounds, &
                          shiftr, shifti )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      integer    ishift, kev, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 bounds(kev+np), ritzr(kev+np), ritzi(kev+np), &
                 shiftr(1), shifti(1)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0, zero = 0.0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy, ssortc
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = mngets
!
!     %----------------------------------------------------%
!     | LM, SM, LR, SR, LI, SI case.                       |
!     | Sort the eigenvalues of H into the desired order   |
!     | and apply the resulting order to BOUNDS.           |
!     | The eigenvalues are sorted so that the wanted part |
!     | are always in the last KEV locations.              |
!     | We first do a pre-processing sort in order to keep |
!     | complex conjugate pairs together                   |
!     %----------------------------------------------------%
!
      if (which == 'LM') then
         call ssortc ('LR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which == 'SM') then
         call ssortc ('SR', .true., kev+np, ritzr, ritzi, bounds)
      else if (which == 'LR') then
         call ssortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which == 'SR') then
         call ssortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which == 'LI') then
         call ssortc ('LM', .true., kev+np, ritzr, ritzi, bounds)
      else if (which == 'SI') then
         call ssortc ('SM', .true., kev+np, ritzr, ritzi, bounds)
      end if
!
      call ssortc (which, .true., kev+np, ritzr, ritzi, bounds)
!
!     %-------------------------------------------------------%
!     | Increase KEV by one if the ( ritzr(np),ritzi(np) )    |
!     | = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) /= zero |
!     | Accordingly decrease NP by one. In other words keep   |
!     | complex conjugate pairs together.                     |
!     %-------------------------------------------------------%
!
      if (       ( ritzr(np+1) - ritzr(np) ) == zero &
           .and. ( ritzi(np+1) + ritzi(np) ) == zero ) then
         np = np - 1
         kev = kev + 1
      end if
!
      if ( ishift == 1 ) then
!
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first        |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when they shifts |
!        | are applied in subroutine snapps.                     |
!        | Be careful and use 'SR' since we want to sort BOUNDS! |
!        %-------------------------------------------------------%
!
         call ssortc ( 'SR', .true., np, bounds, ritzr, ritzi )
      end if
!
      call cpu_time (t1)
      tngets = tngets + (t1 - t0)
!
      if (msglvl > 0) then
         call ivout (logfil, 1, kev, ndigit, '_ngets: KEV is')
         call ivout (logfil, 1, np, ndigit, '_ngets: NP is')
         call svout (logfil, kev+np, ritzr, ndigit, &
              '_ngets: Eigenvalues of current H matrix -- real part')
         call svout (logfil, kev+np, ritzi, ndigit, &
              '_ngets: Eigenvalues of current H matrix -- imag part')
         call svout (logfil, kev+np, bounds, ndigit, &
            '_ngets: Ritz estimates of the current KEV+NP Ritz values')
      end if
!
      return
!
!     %---------------%
!     | End of sngets |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssaitr
!
!\Description:
!  Reverse communication interface for applying NP additional steps to
!  a K step symmetric Arnoldi factorization.
!
!  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
!
!          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
!
!  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
!
!          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
!
!  where OP and B are as in ssaupd.  The B-norm of r_{k+p} is also
!  computed and returned.
!
!\Usage:
!  call ssaitr
!     ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH,
!       IPNTR, WORKD, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!                    This is for the restart phase to force the new
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y,
!                    IPNTR(3) is the pointer into WORK for B * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORK for X,
!                    IPNTR(2) is the pointer into WORK for Y.
!          IDO = 99: done
!          -------------------------------------------------------------
!          When the routine is used in the "shift-and-invert" mode, the
!          vector B * Q is already available and does not need to be
!          recomputed in forming OP * Q.
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of matrix B that defines the
!          semi-inner product for the operator OP.  See ssaupd.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  K       Integer.  (INPUT)
!          Current order of H and the number of columns of V.
!
!  NP      Integer.  (INPUT)
!          Number of additional Arnoldi steps to take.
!
!  MODE    Integer.  (INPUT)
!          Signifies which form for "OP". If MODE=2 then
!          a reduction in the number of B matrix vector multiplies
!          is possible since the B-norm of OP*x is equivalent to
!          the inv(B)-norm of A*x.
!
!  RESID   Real array of length N.  (INPUT/OUTPUT)
!          On INPUT:  RESID contains the residual vector r_{k}.
!          On OUTPUT: RESID contains the residual vector r_{k+p}.
!
!  RNORM   Real scalar.  (INPUT/OUTPUT)
!          On INPUT the B-norm of r_{k}.
!          On OUTPUT the B-norm of the updated residual r_{k+p}.
!
!  V       Real N by K+NP array.  (INPUT/OUTPUT)
!          On INPUT:  V contains the Arnoldi vectors in the first K
!          columns.
!          On OUTPUT: V contains the new NP Arnoldi vectors in the next
!          NP columns.  The first K columns are unchanged.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (K+NP) by 2 array.  (INPUT/OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          with the subdiagonal in the first column starting at H(2,1)
!          and the main diagonal in the second column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORK for
!          vectors used by the Arnoldi iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in the
!                    shift-and-invert mode.  X is the current operand.
!          -------------------------------------------------------------
!
!  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The calling program should not
!          use WORKD as temporary workspace during the iteration !!!!!!
!          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K step Arnoldi factorization. Used to save some
!          computation at the first step.
!          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
!          with the K+NP step Arnoldi factorization.
!
!  INFO    Integer.  (OUTPUT)
!          = 0: Normal exit.
!          > 0: Size of an invariant subspace of OP is found that is
!               less than K + NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     sgetv0  ARPACK routine to generate the initial vector.
!     ivout   ARPACK utility routine that prints integers.
!     smout   ARPACK utility routine that prints matrices.
!     svout   ARPACK utility routine that prints vectors.
!     slamch  LAPACK routine that determines machine constants.
!     slascl  LAPACK routine for careful scaling of a matrix.
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     saxpy   Level 1 BLAS that computes a vector triad.
!     sscal   Level 1 BLAS that scales a vector.
!     scopy   Level 1 BLAS that copies one vector to another .
!     sdot    Level 1 BLAS that computes the scalar product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/93: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2
!
!\Remarks
!  The algorithm implemented is:
!
!  restart = .false.
!  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
!  r_{k} contains the initial residual vector even for k = 0;
!  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
!  computed by the calling program.
!
!  betaj = rnorm ; p_{k+1} = B*r_{k} ;
!  For  j = k+1, ..., k+np  Do
!     1) if ( betaj < tol ) stop or restart depending on j.
!        if ( restart ) generate a new starting vector.
!     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
!        p_{j} = p_{j}/betaj
!     3) r_{j} = OP*v_{j} where OP is defined as in ssaupd
!        For shift-invert mode p_{j} = B*v_{j} is already available.
!        wnorm = || OP*v_{j} ||
!     4) Compute the j-th step residual vector.
!        w_{j} =  V_{j}^T * B * OP * v_{j}
!        r_{j} =  OP*v_{j} - V_{j} * w_{j}
!        alphaj <- j-th component of w_{j}
!        rnorm = || r_{j} ||
!        betaj+1 = rnorm
!        If (rnorm > 0.717*wnorm) accept step and go back to 1)
!     5) Re-orthogonalization step:
!        s = V_{j}'*B*r_{j}
!        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
!        alphaj = alphaj + s_{j};
!     6) Iterative refinement step:
!        If (rnorm1 > 0.717*rnorm) then
!           rnorm = rnorm1
!           accept step and go back to 1)
!        Else
!           rnorm = rnorm1
!           If this is the first time in step 6), go to 5)
!           Else r_{j} lies in the span of V_{j} numerically.
!              Set r_{j} = 0 and rnorm = 0; go to 1)
!        EndIf
!  End Do
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssaitr &
         (ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh, &
          ipntr, workd, info)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, mode, np
      Real &
                 rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Real &
                 h(ldh,2), resid(n), v(ldv,k+np), workd(3*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      logical    first, orth1, orth2, rstart, step3, step4
      integer    i, ierr, ipj, irj, ivj, iter, itry, j, msglvl, &
                 infol, jj
      Real &
                 rnorm1, wnorm, safmin, temp1
      save       orth1, orth2, rstart, step3, step4, &
                 ierr, ipj, irj, ivj, iter, itry, j, msglvl, &
                 rnorm1, safmin, wnorm
!
!     %-----------------------%
!     | Local Array Arguments |
!     %-----------------------%
!
      Real &
                 xtemp(2)
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   saxpy, scopy, sscal, sgemv, sgetv0, svout, smout, &
                 slascl, ivout
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 sdot, snrm2, slamch
      external   sdot, snrm2, slamch
!
!     %-----------------%
!     | Data statements |
!     %-----------------%
!
      data      first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (first) then
         first = .false.
!
!        %--------------------------------%
!        | safmin = safe minimum is such  |
!        | that 1/sfmin does not overflow |
!        %--------------------------------%
!
         safmin = slamch('safmin')
      end if
!
      if (ido == 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call cpu_time (t0)
         msglvl = msaitr
!
!        %------------------------------%
!        | Initial call to this routine |
!        %------------------------------%
!
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
!
!        %--------------------------------%
!        | Pointer to the current step of |
!        | the factorization to build     |
!        %--------------------------------%
!
         j      = k + 1
!
!        %------------------------------------------%
!        | Pointers used for reverse communication  |
!        | when using WORKD.                        |
!        %------------------------------------------%
!
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
!
!     %-------------------------------------------------%
!     | When in reverse communication mode one of:      |
!     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
!     | will be .true.                                  |
!     | STEP3: return from computing OP*v_{j}.          |
!     | STEP4: return from computing B-norm of OP*v_{j} |
!     | ORTH1: return from computing B-norm of r_{j+1}  |
!     | ORTH2: return from computing B-norm of          |
!     |        correction to the residual vector.       |
!     | RSTART: return from OP computations needed by   |
!     |         sgetv0.                                 |
!     %-------------------------------------------------%
!
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
!
!     %------------------------------%
!     | Else this is the first step. |
!     %------------------------------%
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |        A R N O L D I     I T E R A T I O N     L O O P       |
!     |                                                              |
!     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
!     %--------------------------------------------------------------%
!
 1000 continue
!
         if (msglvl > 2) then
            call ivout (logfil, 1, j, ndigit, &
                        '_saitr: generating Arnoldi vector no.')
            call svout (logfil, 1, rnorm, ndigit, &
                        '_saitr: B-norm of the current residual =')
         end if
!
!        %---------------------------------------------------------%
!        | Check for exact zero. Equivalent to determing whether a |
!        | j-step Arnoldi factorization is present.                |
!        %---------------------------------------------------------%
!
         if (rnorm > zero) go to 40
!
!           %---------------------------------------------------%
!           | Invariant subspace found, generate a new starting |
!           | vector which is orthogonal to the current Arnoldi |
!           | basis and continue the iteration.                 |
!           %---------------------------------------------------%
!
            if (msglvl > 0) then
               call ivout (logfil, 1, j, ndigit, &
                           '_saitr: ****** restart at step ******')
            end if
!
!           %---------------------------------------------%
!           | ITRY is the loop variable that controls the |
!           | maximum amount of times that a restart is   |
!           | attempted. NRSTRT is used by stat.h         |
!           %---------------------------------------------%
!
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
!
!           %--------------------------------------%
!           | If in reverse communication mode and |
!           | RSTART = .true. flow returns here.   |
!           %--------------------------------------%
!
            call sgetv0 (ido, bmat, itry, .false., n, j, v, ldv, &
                         resid, rnorm, ipntr, workd, ierr)
            if (ido /= 99) go to 9000
            if (ierr < 0) then
               itry = itry + 1
               if (itry <= 3) go to 20
!
!              %------------------------------------------------%
!              | Give up after several restart attempts.        |
!              | Set INFO to the size of the invariant subspace |
!              | which spans OP and exit.                       |
!              %------------------------------------------------%
!
               info = j - 1
               call cpu_time (t1)
               tsaitr = tsaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
!
   40    continue
!
!        %---------------------------------------------------------%
!        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
!        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
!        | when reciprocating a small RNORM, test against lower    |
!        | machine bound.                                          |
!        %---------------------------------------------------------%
!
         call scopy (n, resid, 1, v(1,j), 1)
         if (rnorm >= safmin) then
             temp1 = one / rnorm
             call sscal (n, temp1, v(1,j), 1)
             call sscal (n, temp1, workd(ipj), 1)
         else
!
!            %-----------------------------------------%
!            | To scale both v_{j} and p_{j} carefully |
!            | use LAPACK routine SLASCL               |
!            %-----------------------------------------%
!
             call slascl ('General', i, i, rnorm, one, n, 1, &
                          v(1,j), n, infol)
             call slascl ('General', i, i, rnorm, one, n, 1, &
                          workd(ipj), n, infol)
         end if
!
!        %------------------------------------------------------%
!        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
!        | Note that this is not quite yet r_{j}. See STEP 4    |
!        %------------------------------------------------------%
!
         step3 = .true.
         nopx  = nopx + 1
         call cpu_time (t2)
         call scopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
!
!        %-----------------------------------%
!        | Exit in order to compute OP*v_{j} |
!        %-----------------------------------%
!
         go to 9000
   50    continue
!
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
!        %-----------------------------------%
!
         call cpu_time (t3)
         tmvopx = tmvopx + (t3 - t2)
!
         step3 = .false.
!
!        %------------------------------------------%
!        | Put another copy of OP*v_{j} into RESID. |
!        %------------------------------------------%
!
         call scopy (n, workd(irj), 1, resid, 1)
!
!        %-------------------------------------------%
!        | STEP 4:  Finish extending the symmetric   |
!        |          Arnoldi to length j. If MODE = 2 |
!        |          then B*OP = B*inv(B)*A = A and   |
!        |          we don't need to compute B*OP.   |
!        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
!        | assumed to have A*v_{j}.                  |
!        %-------------------------------------------%
!
         if (mode == 2) go to 65
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %-------------------------------------%
!           | Exit in order to compute B*OP*v_{j} |
!           %-------------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
              call scopy(n, resid, 1 , workd(ipj), 1)
         end if
   60    continue
!
!        %-----------------------------------%
!        | Back from reverse communication;  |
!        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
!        %-----------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         step4 = .false.
!
!        %-------------------------------------%
!        | The following is needed for STEP 5. |
!        | Compute the B-norm of OP*v_{j}.     |
!        %-------------------------------------%
!
   65    continue
         if (mode == 2) then
!
!           %----------------------------------%
!           | Note that the B-norm of OP*v_{j} |
!           | is the inv(B)-norm of A*v_{j}.   |
!           %----------------------------------%
!
            wnorm = sdot (n, resid, 1, workd(ivj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat == 'G') then
            wnorm = sdot (n, resid, 1, workd(ipj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat == 'I') then
            wnorm = snrm2(n, resid, 1)
         end if
!
!        %-----------------------------------------%
!        | Compute the j-th residual corresponding |
!        | to the j step factorization.            |
!        | Use Classical Gram Schmidt and compute: |
!        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
!        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
!        %-----------------------------------------%
!
!
!        %------------------------------------------%
!        | Compute the j Fourier coefficients w_{j} |
!        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
!        %------------------------------------------%
!
         if (mode /= 2 ) then
            call sgemv('T', n, j, one, v, ldv, workd(ipj), 1, zero, &
                        workd(irj), 1)
         else if (mode == 2) then
            call sgemv('T', n, j, one, v, ldv, workd(ivj), 1, zero, &
                        workd(irj), 1)
         end if
!
!        %--------------------------------------%
!        | Orthgonalize r_{j} against V_{j}.    |
!        | RESID contains OP*v_{j}. See STEP 3. |
!        %--------------------------------------%
!
         call sgemv('N', n, j, -one, v, ldv, workd(irj), 1, one, &
                     resid, 1)
!
!        %--------------------------------------%
!        | Extend H to have j rows and columns. |
!        %--------------------------------------%
!
         h(j,2) = workd(irj + j - 1)
         if (j == 1  .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         end if
         call cpu_time (t4)
!
         orth1 = .true.
         iter  = 0
!
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*r_{j} |
!           %----------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH1 = .true. |
!        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
!        %---------------------------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         orth1 = .false.
!
!        %------------------------------%
!        | Compute the B-norm of r_{j}. |
!        %------------------------------%
!
         if (bmat == 'G') then
            rnorm = sdot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat == 'I') then
            rnorm = snrm2(n, resid, 1)
         end if
!
!        %-----------------------------------------------------------%
!        | STEP 5: Re-orthogonalization / Iterative refinement phase |
!        | Maximum NITER_ITREF tries.                                |
!        |                                                           |
!        |          s      = V_{j}^T * B * r_{j}                     |
!        |          r_{j}  = r_{j} - V_{j}*s                         |
!        |          alphaj = alphaj + s_{j}                          |
!        |                                                           |
!        | The stopping criteria used for iterative refinement is    |
!        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
!        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
!        | Determine if we need to correct the residual. The goal is |
!        | to enforce ||v(:,1:j)^T * r_{j}|| <= eps * || r_{j} ||  |
!        %-----------------------------------------------------------%
!
         if (rnorm > 0.717*wnorm) go to 100
         nrorth = nrorth + 1
!
!        %---------------------------------------------------%
!        | Enter the Iterative refinement phase. If further  |
!        | refinement is necessary, loop back here. The loop |
!        | variable is ITER. Perform a step of Classical     |
!        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
!        %---------------------------------------------------%
!
   80    continue
!
         if (msglvl > 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call svout (logfil, 2, xtemp, ndigit, &
                 '_saitr: re-orthonalization ; wnorm and rnorm are')
         end if
!
!        %----------------------------------------------------%
!        | Compute V_{j}^T * B * r_{j}.                       |
!        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
!        %----------------------------------------------------%
!
         call sgemv ('T', n, j, one, v, ldv, workd(ipj), 1, &
                     zero, workd(irj), 1)
!
!        %----------------------------------------------%
!        | Compute the correction to the residual:      |
!        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
!        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
!        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
!        | H(j,j) is updated.                           |
!        %----------------------------------------------%
!
         call sgemv ('N', n, j, -one, v, ldv, workd(irj), 1, &
                     one, resid, 1)
!
         if (j == 1  .or.  rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workd(irj + j - 1)
!
         orth2 = .true.
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
!
!           %-----------------------------------%
!           | Exit in order to compute B*r_{j}. |
!           | r_{j} is the corrected residual.  |
!           %-----------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
            call scopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
!
!        %---------------------------------------------------%
!        | Back from reverse communication if ORTH2 = .true. |
!        %---------------------------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
!        %-----------------------------------------------------%
!        | Compute the B-norm of the corrected residual r_{j}. |
!        %-----------------------------------------------------%
!
         if (bmat == 'G') then
             rnorm1 = sdot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat == 'I') then
             rnorm1 = snrm2(n, resid, 1)
         end if
!
         if (msglvl > 0 .and. iter .gt. 0) then
            call ivout (logfil, 1, j, ndigit, &
                 '_saitr: Iterative refinement for Arnoldi residual')
            if (msglvl > 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call svout (logfil, 2, xtemp, ndigit, &
                 '_saitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
!
!        %-----------------------------------------%
!        | Determine if we need to perform another |
!        | step of re-orthogonalization.           |
!        %-----------------------------------------%
!
         if (rnorm1 > 0.717*rnorm) then
!
!           %--------------------------------%
!           | No need for further refinement |
!           %--------------------------------%
!
            rnorm = rnorm1
!
         else
!
!           %-------------------------------------------%
!           | Another step of iterative refinement step |
!           | is required. NITREF is used by stat.h     |
!           %-------------------------------------------%
!
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter <= 1) go to 80
!
!           %-------------------------------------------------%
!           | Otherwise RESID is numerically in the span of V |
!           %-------------------------------------------------%
!
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
!
!        %----------------------------------------------%
!        | Branch here directly if iterative refinement |
!        | wasn't necessary or after at most NITER_REF  |
!        | steps of iterative refinement.               |
!        %----------------------------------------------%
!
  100    continue
!
         rstart = .false.
         orth2  = .false.
!
         call cpu_time (t5)
         titref = titref + (t5 - t4)
!
!        %----------------------------------------------------------%
!        | Make sure the last off-diagonal element is non negative  |
!        | If not perform a similarity transformation on H(1:j,1:j) |
!        | and scale v(:,j) by -1.                                  |
!        %----------------------------------------------------------%
!
         if (h(j,1) < zero) then
            h(j,1) = -h(j,1)
            if ( j < k+np) then
               call sscal(n, -one, v(1,j+1), 1)
            else
               call sscal(n, -one, resid, 1)
            end if
         end if
!
!        %------------------------------------%
!        | STEP 6: Update  j = j+1;  Continue |
!        %------------------------------------%
!
         j = j + 1
         if (j > k+np) then
            call cpu_time (t1)
            tsaitr = tsaitr + (t1 - t0)
            ido = 99
!
            if (msglvl > 1) then
               call svout (logfil, k+np, h(1,2), ndigit, &
               '_saitr: main diagonal of matrix H of step K+NP.')
               if (k+np > 1) then
               call svout (logfil, k+np-1, h(2,1), ndigit, &
               '_saitr: sub diagonal of matrix H of step K+NP.')
               end if
            end if
!
            go to 9000
         end if
!
!        %--------------------------------------------------------%
!        | Loop back to extend the factorization by another step. |
!        %--------------------------------------------------------%
!
      go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 9000 continue
      return
!
!     %---------------%
!     | End of ssaitr |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssapps
!
!\Description:
!  Given the Arnoldi factorization
!
!     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,
!
!  apply NP shifts implicitly resulting in
!
!     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q
!
!  where Q is an orthogonal matrix of order KEV+NP. Q is the product of
!  rotations resulting from the NP bulge chasing sweeps.  The updated Arnoldi
!  factorization becomes:
!
!     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.
!
!\Usage:
!  call ssapps
!     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ, WORKD )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Problem size, i.e. dimension of matrix A.
!
!  KEV     Integer.  (INPUT)
!          INPUT: KEV+NP is the size of the input matrix H.
!          OUTPUT: KEV is the size of the updated matrix HNEW.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be applied.
!
!  SHIFT   Real array of length NP.  (INPUT)
!          The shifts to be applied.
!
!  V       Real N by (KEV+NP) array.  (INPUT/OUTPUT)
!          INPUT: V contains the current KEV+NP Arnoldi vectors.
!          OUTPUT: VNEW = V(1:n,1:KEV); the updated Arnoldi vectors
!          are in the first KEV columns of V.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (KEV+NP) by 2 array.  (INPUT/OUTPUT)
!          INPUT: H contains the symmetric tridiagonal matrix of the
!          Arnoldi factorization with the subdiagonal in the 1st column
!          starting at H(2,1) and the main diagonal in the 2nd column.
!          OUTPUT: H contains the updated tridiagonal matrix in the
!          KEV leading submatrix.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RESID   Real array of length (N).  (INPUT/OUTPUT)
!          INPUT: RESID contains the the residual vector r_{k+p}.
!          OUTPUT: RESID is the updated residual vector rnew_{k}.
!
!  Q       Real KEV+NP by KEV+NP work array.  (WORKSPACE)
!          Work array used to accumulate the rotations during the bulge
!          chase sweep.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKD   Real work array of length 2*N.  (WORKSPACE)
!          Distributed array used in the application of the accumulated
!          orthogonal matrix Q.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!
!\Routines called:
!     ivout   ARPACK utility routine that prints integers.
!     svout   ARPACK utility routine that prints vectors.
!     slamch  LAPACK routine that determines machine constants.
!     slartg  LAPACK Givens rotation construction routine.
!     slacpy  LAPACK matrix copy routine.
!     slaset  LAPACK matrix initialization routine.
!     sgemv   Level 2 BLAS routine for matrix vector multiplication.
!     saxpy   Level 1 BLAS that computes a vector triad.
!     scopy   Level 1 BLAS that copies one vector to another.
!     sscal   Level 1 BLAS that scales a vector.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/16/93: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: sapps.F   SID: 2.6   DATE OF SID: 3/28/97   RELEASE: 2
!
!\Remarks
!  1. In this version, each shift is applied to all the subblocks of
!     the tridiagonal matrix H and not just to the submatrix that it
!     comes from. This routine assumes that the subdiagonal elements
!     of H that are stored in h(1:kev+np,1) are nonegative upon input
!     and enforce this condition upon output. This version incorporates
!     deflation. See code for documentation.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssapps &
         ( n, kev, np, shift, v, ldv, h, ldh, resid, q, ldq, workd )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    kev, ldh, ldq, ldv, n, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 h(ldh,2), q(ldq,kev+np), resid(n), shift(np), &
                 v(ldv,kev+np), workd(2*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, iend, istart, itop, j, jj, kplusp, msglvl
      logical    first
      Real &
                 a1, a2, a3, a4, big, c, epsmch, f, g, r, s
      save       epsmch, first
!
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   saxpy, scopy, sscal, slacpy, slartg, slaset, svout, &
                 ivout, sgemv
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slamch
      external   slamch
!
!     %----------------------%
!     | Intrinsics Functions |
!     %----------------------%
!
      intrinsic  abs
!
!     %----------------%
!     | Data statments |
!     %----------------%
!
      data       first / .true. /
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (first) then
         epsmch = slamch('Epsilon-Machine')
         first = .false.
      end if
      itop = 1
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = msapps
!
      kplusp = kev + np
!
!     %----------------------------------------------%
!     | Initialize Q to the identity matrix of order |
!     | kplusp used to accumulate the rotations.     |
!     %----------------------------------------------%
!
      call slaset ('All', kplusp, kplusp, zero, one, q, ldq)
!
!     %----------------------------------------------%
!     | Quick return if there are no shifts to apply |
!     %----------------------------------------------%
!
      if (np == 0) go to 9000
!
!     %----------------------------------------------------------%
!     | Apply the np shifts implicitly. Apply each shift to the  |
!     | whole matrix and not just to the submatrix from which it |
!     | comes.                                                   |
!     %----------------------------------------------------------%
!
      do 90 jj = 1, np
!
         istart = itop
!
!        %----------------------------------------------------------%
!        | Check for splitting and deflation. Currently we consider |
!        | an off-diagonal element h(i+1,1) negligible if           |
!        |         h(i+1,1) <= epsmch*( |h(i,2)| + |h(i+1,2)| )   |
!        | for i=1:KEV+NP-1.                                        |
!        | If above condition tests true then we set h(i+1,1) = 0.  |
!        | Note that h(1:KEV+NP,1) are assumed to be non negative.  |
!        %----------------------------------------------------------%
!
   20    continue
!
!        %------------------------------------------------%
!        | The following loop exits early if we encounter |
!        | a negligible off diagonal element.             |
!        %------------------------------------------------%
!
         do 30 i = istart, kplusp-1
            big   = abs(h(i,2)) + abs(h(i+1,2))
            if (h(i+1,1) <= epsmch*big) then
               if (msglvl > 0) then
                  call ivout (logfil, 1, i, ndigit, &
                       '_sapps: deflation at row/column no.')
                  call ivout (logfil, 1, jj, ndigit, &
                       '_sapps: occured before shift number.')
                  call svout (logfil, 1, h(i+1,1), ndigit, &
                       '_sapps: the corresponding off diagonal element')
               end if
               h(i+1,1) = zero
               iend = i
               go to 40
            end if
   30    continue
         iend = kplusp
   40    continue
!
         if (istart < iend) then
!
!           %--------------------------------------------------------%
!           | Construct the plane rotation G'(istart,istart+1,theta) |
!           | that attempts to drive h(istart+1,1) to zero.          |
!           %--------------------------------------------------------%
!
             f = h(istart,2) - shift(jj)
             g = h(istart+1,1)
             call slartg (f, g, c, s, r)
!
!            %-------------------------------------------------------%
!            | Apply rotation to the left and right of H;            |
!            | H <- G' * H * G,  where G = G(istart,istart+1,theta). |
!            | This will create a "bulge".                           |
!            %-------------------------------------------------------%
!
             a1 = c*h(istart,2)   + s*h(istart+1,1)
             a2 = c*h(istart+1,1) + s*h(istart+1,2)
             a4 = c*h(istart+1,2) - s*h(istart+1,1)
             a3 = c*h(istart+1,1) - s*h(istart,2)
             h(istart,2)   = c*a1 + s*a2
             h(istart+1,2) = c*a4 - s*a3
             h(istart+1,1) = c*a3 + s*a4
!
!            %----------------------------------------------------%
!            | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!            %----------------------------------------------------%
!
             do 60 j = 1, min(istart+jj,kplusp)
                a1            =   c*q(j,istart) + s*q(j,istart+1)
                q(j,istart+1) = - s*q(j,istart) + c*q(j,istart+1)
                q(j,istart)   = a1
   60        continue
!
!
!            %----------------------------------------------%
!            | The following loop chases the bulge created. |
!            | Note that the previous rotation may also be  |
!            | done within the following loop. But it is    |
!            | kept separate to make the distinction among  |
!            | the bulge chasing sweeps and the first plane |
!            | rotation designed to drive h(istart+1,1) to  |
!            | zero.                                        |
!            %----------------------------------------------%
!
             do 70 i = istart+1, iend-1
!
!               %----------------------------------------------%
!               | Construct the plane rotation G'(i,i+1,theta) |
!               | that zeros the i-th bulge that was created   |
!               | by G(i-1,i,theta). g represents the bulge.   |
!               %----------------------------------------------%
!
                f = h(i,1)
                g = s*h(i+1,1)
!
!               %----------------------------------%
!               | Final update with G(i-1,i,theta) |
!               %----------------------------------%
!
                h(i+1,1) = c*h(i+1,1)
                call slartg (f, g, c, s, r)
!
!               %-------------------------------------------%
!               | The following ensures that h(1:iend-1,1), |
!               | the first iend-2 off diagonal of elements |
!               | H, remain non negative.                   |
!               %-------------------------------------------%
!
                if (r < zero) then
                   r = -r
                   c = -c
                   s = -s
                end if
!
!               %--------------------------------------------%
!               | Apply rotation to the left and right of H; |
!               | H <- G * H * G',  where G = G(i,i+1,theta) |
!               %--------------------------------------------%
!
                h(i,1) = r
!
                a1 = c*h(i,2)   + s*h(i+1,1)
                a2 = c*h(i+1,1) + s*h(i+1,2)
                a3 = c*h(i+1,1) - s*h(i,2)
                a4 = c*h(i+1,2) - s*h(i+1,1)
!
                h(i,2)   = c*a1 + s*a2
                h(i+1,2) = c*a4 - s*a3
                h(i+1,1) = c*a3 + s*a4
!
!               %----------------------------------------------------%
!               | Accumulate the rotation in the matrix Q;  Q <- Q*G |
!               %----------------------------------------------------%
!
                do 50 j = 1, min( i+jj, kplusp )
                   a1       =   c*q(j,i) + s*q(j,i+1)
                   q(j,i+1) = - s*q(j,i) + c*q(j,i+1)
                   q(j,i)   = a1
   50           continue
!
   70        continue
!
         end if
!
!        %--------------------------%
!        | Update the block pointer |
!        %--------------------------%
!
         istart = iend + 1
!
!        %------------------------------------------%
!        | Make sure that h(iend,1) is non-negative |
!        | If not then set h(iend,1) <-- -h(iend,1) |
!        | and negate the last column of Q.         |
!        | We have effectively carried out a        |
!        | similarity on transformation H           |
!        %------------------------------------------%
!
         if (h(iend,1) < zero) then
             h(iend,1) = -h(iend,1)
             call sscal(kplusp, -one, q(1,iend), 1)
         end if
!
!        %--------------------------------------------------------%
!        | Apply the same shift to the next block if there is any |
!        %--------------------------------------------------------%
!
         if (iend < kplusp) go to 20
!
!        %-----------------------------------------------------%
!        | Check if we can increase the the start of the block |
!        %-----------------------------------------------------%
!
         do 80 i = itop, kplusp-1
            if (h(i+1,1) > zero) go to 90
            itop  = itop + 1
   80    continue
!
!        %-----------------------------------%
!        | Finished applying the jj-th shift |
!        %-----------------------------------%
!
   90 continue
!
!     %------------------------------------------%
!     | All shifts have been applied. Check for  |
!     | more possible deflation that might occur |
!     | after the last shift is applied.         |
!     %------------------------------------------%
!
      do i = itop, kplusp-1
         big   = abs(h(i,2)) + abs(h(i+1,2))
         if (h(i+1,1) <= epsmch*big) then
            if (msglvl > 0) then
               call ivout (logfil, 1, i, ndigit, &
                    '_sapps: deflation at row/column no.')
               call svout (logfil, 1, h(i+1,1), ndigit, &
                    '_sapps: the corresponding off diagonal element')
            end if
            h(i+1,1) = zero
         end if
      end do
!
!     %-------------------------------------------------%
!     | Compute the (kev+1)-st column of (V*Q) and      |
!     | temporarily store the result in WORKD(N+1:2*N). |
!     | This is not necessary if h(kev+1,1) = 0.         |
!     %-------------------------------------------------%
!
      if ( h(kev+1,1) > zero ) &
         call sgemv ('N', n, kplusp, one, v, ldv, &
                      q(1,kev+1), 1, zero, workd(n+1), 1)
!
!     %-------------------------------------------------------%
!     | Compute column 1 to kev of (V*Q) in backward order    |
!     | taking advantage that Q is an upper triangular matrix |
!     | with lower bandwidth np.                              |
!     | Place results in v(:,kplusp-kev:kplusp) temporarily.  |
!     %-------------------------------------------------------%
!
      do 130 i = 1, kev
         call sgemv ('N', n, kplusp-i+1, one, v, ldv, &
                     q(1,kev-i+1), 1, zero, workd, 1)
         call scopy (n, workd, 1, v(1,kplusp-i+1), 1)
  130 continue
!
!     %-------------------------------------------------%
!     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
!     %-------------------------------------------------%
!
      call slacpy ('All', n, kev, v(1,np+1), ldv, v, ldv)
!
!     %--------------------------------------------%
!     | Copy the (kev+1)-st column of (V*Q) in the |
!     | appropriate place if h(kev+1,1) /= zero. |
!     %--------------------------------------------%
!
      if ( h(kev+1,1) > zero ) &
           call scopy (n, workd(n+1), 1, v(1,kev+1), 1)
!
!     %-------------------------------------%
!     | Update the residual vector:         |
!     |    r <- sigmak*r + betak*v(:,kev+1) |
!     | where                               |
!     |    sigmak = (e_{kev+p}'*Q)*e_{kev}  |
!     |    betak = e_{kev+1}'*H*e_{kev}     |
!     %-------------------------------------%
!
      call sscal (n, q(kplusp,kev), resid, 1)
      if (h(kev+1,1) > zero) &
         call saxpy (n, h(kev+1,1), v(1,kev+1), 1, resid, 1)
!
      if (msglvl > 1) then
         call svout (logfil, 1, q(kplusp,kev), ndigit, &
            '_sapps: sigmak of the updated residual vector')
         call svout (logfil, 1, h(kev+1,1), ndigit, &
            '_sapps: betak of the updated residual vector')
         call svout (logfil, kev, h(1,2), ndigit, &
            '_sapps: updated main diagonal of H for next iteration')
         if (kev > 1) then
         call svout (logfil, kev-1, h(2,1), ndigit, &
            '_sapps: updated sub diagonal of H for next iteration')
         end if
      end if
!
      call cpu_time (t1)
      tsapps = tsapps + (t1 - t0)
!
 9000 continue
      return
!
!     %---------------%
!     | End of ssapps |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssaup2
!
!\Description:
!  Intermediate level interface called by ssaupd.
!
!\Usage:
!  call ssaup2
!     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
!       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL,
!       IPNTR, WORKD, INFO )
!
!\Arguments
!
!  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in ssaupd.
!  MODE, ISHIFT, MXITER: see the definition of IPARAM in ssaupd.
!
!  NP      Integer.  (INPUT/OUTPUT)
!          Contains the number of implicit shifts to apply during
!          each Arnoldi/Lanczos iteration.
!          If ISHIFT=1, NP is adjusted dynamically at each iteration
!          to accelerate convergence and prevent stagnation.
!          This is also roughly equal to the number of matrix-vector
!          products (involving the operator OP) per Arnoldi iteration.
!          The logic for adjusting is contained within the current
!          subroutine.
!          If ISHIFT=0, NP is the number of shifts the user needs
!          to provide via reverse comunication. 0 < NP < NCV-NEV.
!          NP may be less than NCV-NEV since a leading block of the current
!          upper Tridiagonal matrix has split off and contains "unwanted"
!          Ritz values.
!          Upon termination of the IRA iteration, NP contains the number
!          of "converged" wanted Ritz values.
!
!  IUPD    Integer.  (INPUT)
!          IUPD .EQ. 0: use explicit restart instead implicit update.
!          IUPD .NE. 0: use implicit update.
!
!  V       Real N by (NEV+NP) array.  (INPUT/OUTPUT)
!          The Lanczos basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  H       Real (NEV+NP) by 2 array.  (OUTPUT)
!          H is used to store the generated symmetric tridiagonal matrix
!          The subdiagonal is stored in the first column of H starting
!          at H(2,1).  The main diagonal is stored in the second column
!          of H starting at H(1,2). If ssaup2 converges store the
!          B-norm of the final residual vector in H(1,1).
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  RITZ    Real array of length NEV+NP.  (OUTPUT)
!          RITZ(1:NEV) contains the computed Ritz values of OP.
!
!  BOUNDS  Real array of length NEV+NP.  (OUTPUT)
!          BOUNDS(1:NEV) contain the error bounds corresponding to RITZ.
!
!  Q       Real (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
!          Private (replicated) work array used to accumulate the
!          rotation in the shift application step.
!
!  LDQ     Integer.  (INPUT)
!          Leading dimension of Q exactly as declared in the calling
!          program.
!
!  WORKL   Real array of length at least 3*(NEV+NP).  (INPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  It is used in the computation of the
!          tridiagonal eigenvalue problem, the calculation and
!          application of the shifts and convergence checking.
!          If ISHIFT .EQ. O and IDO .EQ. 3, the first NP locations
!          of WORKL are used in reverse communication to hold the user
!          supplied shifts.
!
!  IPNTR   Integer array of length 3.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD for
!          vectors used by the Lanczos iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X.
!          IPNTR(2): pointer to the current result vector Y.
!          IPNTR(3): pointer to the vector B * X when used in one of
!                    the spectral transformation modes.  X is the current
!                    operand.
!          -------------------------------------------------------------
!
!  WORKD   Real work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Lanczos iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration !!!!!!!!!!
!          See Data Distribution Note in ssaupd.
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =     0: Normal return.
!          =     1: All possible eigenvalues of OP has been found.
!                   NP returns the size of the invariant subspace
!                   spanning the operator OP.
!          =     2: No shifts could be applied.
!          =    -8: Error return from trid. eigenvalue calculation;
!                   This should never happen.
!          =    -9: Starting vector is zero.
!          = -9999: Could not build an Lanczos factorization.
!                   Size that was built in returned in NP.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Routines called:
!     sgetv0  ARPACK initial vector generation routine.
!     ssaitr  ARPACK Lanczos factorization routine.
!     ssapps  ARPACK application of implicit shifts routine.
!     ssconv  ARPACK convergence of Ritz values routine.
!     sseigt  ARPACK compute Ritz values and error bounds routine.
!     ssgets  ARPACK reorder Ritz values and error bounds routine.
!     ssortr  ARPACK sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     svout   ARPACK utility routine that prints vectors.
!     slamch  LAPACK routine that determines machine constants.
!     scopy   Level 1 BLAS that copies one vector to another.
!     sdot    Level 1 BLAS that computes the scalar product of two vectors.
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sscal   Level 1 BLAS that scales a vector.
!     sswap   Level 1 BLAS that swaps two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/15/93: Version ' 2.4'
!     xx/xx/95: Version ' 2.4'.  (R.B. Lehoucq)
!
!\SCCS Information: @(#)
! FILE: saup2.F   SID: 2.7   DATE OF SID: 5/19/98   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssaup2 &
         ( ido, bmat, n, which, nev, np, tol, resid, mode, iupd, &
           ishift, mxiter, v, ldv, h, ldh, ritz, bounds, &
           q, ldq, workl, ipntr, workd, info )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1, which*2
      integer    ido, info, ishift, iupd, ldh, ldq, ldv, mxiter, &
                 n, mode, nev, np
      Real &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    ipntr(3)
      Real &
                 bounds(nev+np), h(ldh,2), q(ldq,nev+np), resid(n), &
                 ritz(nev+np), v(ldv,nev+np), workd(3*n), &
                 workl(3*(nev+np))
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  wprime*2
      logical    cnorm, getv0, initv, update, ushift
      integer    ierr, iter, j, kplusp, msglvl, nconv, nevbef, nev0, &
                 np0, nptemp, nevd2, nevm2, kp(3)
      Real &
                 rnorm, temp, eps23
      save       cnorm, getv0, initv, update, ushift, &
                 iter, kplusp, msglvl, nconv, nev0, np0, &
                 rnorm, eps23
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy, sgetv0, ssaitr, sscal, ssconv, sseigt, ssgets, &
                 ssapps, ssortr, svout, ivout, sswap
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 sdot, snrm2, slamch
      external   sdot, snrm2, slamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido == 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call cpu_time (t0)
         msglvl = msaup2
!
!        %---------------------------------%
!        | Set machine dependent constant. |
!        %---------------------------------%
!
         eps23 = slamch('Epsilon-Machine')
         eps23 = eps23**(2.0E+0/3.0E+0)
!
!        %-------------------------------------%
!        | nev0 and np0 are integer variables  |
!        | hold the initial values of NEV & NP |
!        %-------------------------------------%
!
         nev0   = nev
         np0    = np
!
!        %-------------------------------------%
!        | kplusp is the bound on the largest  |
!        |        Lanczos factorization built. |
!        | nconv is the current number of      |
!        |        "converged" eigenvlues.      |
!        | iter is the counter on the current  |
!        |      iteration step.                |
!        %-------------------------------------%
!
         kplusp = nev0 + np0
         nconv  = 0
         iter   = 0
!
!        %--------------------------------------------%
!        | Set flags for computing the first NEV steps |
!        | of the Lanczos factorization.              |
!        %--------------------------------------------%
!
         getv0    = .true.
         update   = .false.
         ushift   = .false.
         cnorm    = .false.
!
         if (info /= 0) then
!
!        %--------------------------------------------%
!        | User provides the initial residual vector. |
!        %--------------------------------------------%
!
            initv = .true.
            info  = 0
         else
            initv = .false.
         end if
      end if
!
!     %---------------------------------------------%
!     | Get a possibly random starting vector and   |
!     | force it into the range of the operator OP. |
!     %---------------------------------------------%
!
   10 continue
!
      if (getv0) then
         call sgetv0 (ido, bmat, 1, initv, n, 1, v, ldv, resid, rnorm, &
                      ipntr, workd, info)
!
         if (ido /= 99) go to 9000
!
         if (rnorm == zero) then
!
!           %-----------------------------------------%
!           | The initial vector is zero. Error exit. |
!           %-----------------------------------------%
!
            info = -9
            go to 1200
         end if
         getv0 = .false.
         ido  = 0
      end if
!
!     %------------------------------------------------------------%
!     | Back from reverse communication: continue with update step |
!     %------------------------------------------------------------%
!
      if (update) go to 20
!
!     %-------------------------------------------%
!     | Back from computing user specified shifts |
!     %-------------------------------------------%
!
      if (ushift) go to 50
!
!     %-------------------------------------%
!     | Back from computing residual norm   |
!     | at the end of the current iteration |
!     %-------------------------------------%
!
      if (cnorm)  go to 100
!
!     %----------------------------------------------------------%
!     | Compute the first NEV steps of the Lanczos factorization |
!     %----------------------------------------------------------%
!
      call ssaitr (ido, bmat, n, 0, nev0, mode, resid, rnorm, v, ldv, &
                   h, ldh, ipntr, workd, info)
!
!     %---------------------------------------------------%
!     | ido /= 99 implies use of reverse communication  |
!     | to compute operations involving OP and possibly B |
!     %---------------------------------------------------%
!
      if (ido /= 99) go to 9000
!
      if (info > 0) then
!
!        %-----------------------------------------------------%
!        | ssaitr was unable to build an Lanczos factorization |
!        | of length NEV0. INFO is returned with the size of   |
!        | the factorization built. Exit main loop.            |
!        %-----------------------------------------------------%
!
         np   = info
         mxiter = iter
         info = -9999
         go to 1200
      end if
!
!     %--------------------------------------------------------------%
!     |                                                              |
!     |           M A I N  LANCZOS  I T E R A T I O N  L O O P       |
!     |           Each iteration implicitly restarts the Lanczos     |
!     |           factorization in place.                            |
!     |                                                              |
!     %--------------------------------------------------------------%
!
 1000 continue
!
         iter = iter + 1
!
         if (msglvl > 0) then
            call ivout (logfil, 1, iter, ndigit, &
                 '_saup2: **** Start of major iteration number ****')
         end if
         if (msglvl > 1) then
            call ivout (logfil, 1, nev, ndigit, &
           '_saup2: The length of the current Lanczos factorization')
            call ivout (logfil, 1, np, ndigit, &
                 '_saup2: Extend the Lanczos factorization by')
         end if
!
!        %------------------------------------------------------------%
!        | Compute NP additional steps of the Lanczos factorization. |
!        %------------------------------------------------------------%
!
         ido = 0
   20    continue
         update = .true.
!
         call ssaitr (ido, bmat, n, nev, np, mode, resid, rnorm, v, &
                      ldv, h, ldh, ipntr, workd, info)
!
!        %---------------------------------------------------%
!        | ido /= 99 implies use of reverse communication  |
!        | to compute operations involving OP and possibly B |
!        %---------------------------------------------------%
!
         if (ido /= 99) go to 9000
!
         if (info > 0) then
!
!           %-----------------------------------------------------%
!           | ssaitr was unable to build an Lanczos factorization |
!           | of length NEV0+NP0. INFO is returned with the size  |
!           | of the factorization built. Exit main loop.         |
!           %-----------------------------------------------------%
!
            np = info
            mxiter = iter
            info = -9999
            go to 1200
         end if
         update = .false.

         if (msglvl > 1) then
            call svout (logfil, 1, rnorm, ndigit, &
                 '_saup2: Current B-norm of residual for factorization')
         end if
!
!        %--------------------------------------------------------%
!        | Compute the eigenvalues and corresponding error bounds |
!        | of the current symmetric tridiagonal matrix.           |
!        %--------------------------------------------------------%
!
         call sseigt (rnorm, kplusp, h, ldh, ritz, bounds, workl, ierr)
!
         if (ierr /= 0) then
            info = -8
            go to 1200
         end if
!
!        %----------------------------------------------------%
!        | Make a copy of eigenvalues and corresponding error |
!        | bounds obtained from _seigt.                       |
!        %----------------------------------------------------%
!
         call scopy(kplusp, ritz, 1, workl(kplusp+1), 1)
         call scopy(kplusp, bounds, 1, workl(2*kplusp+1), 1)
!
!        %---------------------------------------------------%
!        | Select the wanted Ritz values and their bounds    |
!        | to be used in the convergence test.               |
!        | The selection is based on the requested number of |
!        | eigenvalues instead of the current NEV and NP to  |
!        | prevent possible misconvergence.                  |
!        | * Wanted Ritz values := RITZ(NP+1:NEV+NP)         |
!        | * Shifts := RITZ(1:NP) := WORKL(1:NP)             |
!        %---------------------------------------------------%
!
         nev = nev0
         np = np0
         call ssgets (ishift, which, nev, np, ritz, bounds, workl)
!
!        %-------------------%
!        | Convergence test. |
!        %-------------------%
!
         call scopy (nev, bounds(np+1), 1, workl(np+1), 1)
         call ssconv (nev, ritz(np+1), workl(np+1), tol, nconv)
!
         if (msglvl > 2) then
            kp(1) = nev
            kp(2) = np
            kp(3) = nconv
            call ivout (logfil, 3, kp, ndigit, &
                        '_saup2: NEV, NP, NCONV are')
            call svout (logfil, kplusp, ritz, ndigit, &
                 '_saup2: The eigenvalues of H')
            call svout (logfil, kplusp, bounds, ndigit, &
                '_saup2: Ritz estimates of the current NCV Ritz values')
         end if
!
!        %---------------------------------------------------------%
!        | Count the number of unwanted Ritz values that have zero |
!        | Ritz estimates. If any Ritz estimates are equal to zero |
!        | then a leading block of H of order equal to at least    |
!        | the number of Ritz values with zero Ritz estimates has  |
!        | split off. None of these Ritz values may be removed by  |
!        | shifting. Decrease NP the number of shifts to apply. If |
!        | no shifts may be applied, then prepare to exit          |
!        %---------------------------------------------------------%
!
         nptemp = np
         do 30 j=1, nptemp
            if (bounds(j) == zero) then
               np = np - 1
               nev = nev + 1
            end if
 30      continue
!
         if ( (nconv >= nev0) .or. &
              (iter > mxiter) .or. &
              (np == 0) ) then
!
!           %------------------------------------------------%
!           | Prepare to exit. Put the converged Ritz values |
!           | and corresponding bounds in RITZ(1:NCONV) and  |
!           | BOUNDS(1:NCONV) respectively. Then sort. Be    |
!           | careful when NCONV > NP since we don't want to |
!           | swap overlapping locations.                    |
!           %------------------------------------------------%
!
            if (which == 'BE') then
!
!              %-----------------------------------------------------%
!              | Both ends of the spectrum are requested.            |
!              | Sort the eigenvalues into algebraically decreasing  |
!              | order first then swap low end of the spectrum next  |
!              | to high end in appropriate locations.               |
!              | NOTE: when np < floor(nev/2) be careful not to swap |
!              | overlapping locations.                              |
!              %-----------------------------------------------------%
!
               wprime = 'SA'
               call ssortr (wprime, .true., kplusp, ritz, bounds)
               nevd2 = nev0 / 2
               nevm2 = nev0 - nevd2
               if ( nev > 1 ) then
                  call sswap ( min(nevd2,np), ritz(nevm2+1), 1, &
                       ritz( max(kplusp-nevd2+1,kplusp-np+1) ), 1)
                  call sswap ( min(nevd2,np), bounds(nevm2+1), 1, &
                       bounds( max(kplusp-nevd2+1,kplusp-np+1)), 1)
               end if
!
            else
!
!              %--------------------------------------------------%
!              | LM, SM, LA, SA case.                             |
!              | Sort the eigenvalues of H into the an order that |
!              | is opposite to WHICH, and apply the resulting    |
!              | order to BOUNDS.  The eigenvalues are sorted so  |
!              | that the wanted part are always within the first |
!              | NEV locations.                                   |
!              %--------------------------------------------------%
!
               if (which == 'LM') wprime = 'SM'
               if (which == 'SM') wprime = 'LM'
               if (which == 'LA') wprime = 'SA'
               if (which == 'SA') wprime = 'LA'
!
               call ssortr (wprime, .true., kplusp, ritz, bounds)
!
            end if
!
!           %--------------------------------------------------%
!           | Scale the Ritz estimate of each Ritz value       |
!           | by 1 / max(eps23,magnitude of the Ritz value).   |
!           %--------------------------------------------------%
!
            do 35 j = 1, nev0
               temp = max( eps23, abs(ritz(j)) )
               bounds(j) = bounds(j)/temp
 35         continue
!
!           %----------------------------------------------------%
!           | Sort the Ritz values according to the scaled Ritz  |
!           | esitmates.  This will push all the converged ones  |
!           | towards the front of ritzr, ritzi, bounds          |
!           | (in the case when NCONV < NEV.)                    |
!           %----------------------------------------------------%
!
            wprime = 'LA'
            call ssortr(wprime, .true., nev0, bounds, ritz)
!
!           %----------------------------------------------%
!           | Scale the Ritz estimate back to its original |
!           | value.                                       |
!           %----------------------------------------------%
!
            do 40 j = 1, nev0
                temp = max( eps23, abs(ritz(j)) )
                bounds(j) = bounds(j)*temp
 40         continue
!
!           %--------------------------------------------------%
!           | Sort the "converged" Ritz values again so that   |
!           | the "threshold" values and their associated Ritz |
!           | estimates appear at the appropriate position in  |
!           | ritz and bound.                                  |
!           %--------------------------------------------------%
!
            if (which == 'BE') then
!
!              %------------------------------------------------%
!              | Sort the "converged" Ritz values in increasing |
!              | order.  The "threshold" values are in the      |
!              | middle.                                        |
!              %------------------------------------------------%
!
               wprime = 'LA'
               call ssortr(wprime, .true., nconv, ritz, bounds)
!
            else
!
!              %----------------------------------------------%
!              | In LM, SM, LA, SA case, sort the "converged" |
!              | Ritz values according to WHICH so that the   |
!              | "threshold" value appears at the front of    |
!              | ritz.                                        |
!              %----------------------------------------------%

               call ssortr(which, .true., nconv, ritz, bounds)
!
            end if
!
!           %------------------------------------------%
!           |  Use h( 1,1 ) as storage to communicate  |
!           |  rnorm to _seupd if needed               |
!           %------------------------------------------%
!
            h(1,1) = rnorm
!
            if (msglvl > 1) then
               call svout (logfil, kplusp, ritz, ndigit, &
                  '_saup2: Sorted Ritz values.')
               call svout (logfil, kplusp, bounds, ndigit, &
                  '_saup2: Sorted ritz estimates.')
            end if
!
!           %------------------------------------%
!           | Max iterations have been exceeded. |
!           %------------------------------------%
!
            if (iter > mxiter .and. nconv < nev) info = 1
!
!           %---------------------%
!           | No shifts to apply. |
!           %---------------------%
!
            if (np == 0 .and. nconv < nev0) info = 2
!
            np = nconv
            go to 1100
!
         else if (nconv < nev .and. ishift == 1) then
!
!           %---------------------------------------------------%
!           | Do not have all the requested eigenvalues yet.    |
!           | To prevent possible stagnation, adjust the number |
!           | of Ritz values and the shifts.                    |
!           %---------------------------------------------------%
!
            nevbef = nev
            nev = nev + min (nconv, np/2)
            if (nev == 1 .and. kplusp >= 6) then
               nev = kplusp / 2
            else if (nev == 1 .and. kplusp > 2) then
               nev = 2
            end if
            np  = kplusp - nev
!
!           %---------------------------------------%
!           | If the size of NEV was just increased |
!           | resort the eigenvalues.               |
!           %---------------------------------------%
!
            if (nevbef < nev) &
               call ssgets (ishift, which, nev, np, ritz, bounds, &
                    workl)
!
         end if
!
         if (msglvl > 0) then
            call ivout (logfil, 1, nconv, ndigit, &
                 '_saup2: no. of "converged" Ritz values at this iter.')
            if (msglvl > 1) then
               kp(1) = nev
               kp(2) = np
               call ivout (logfil, 2, kp, ndigit, &
                    '_saup2: NEV and NP are')
               call svout (logfil, nev, ritz(np+1), ndigit, &
                    '_saup2: "wanted" Ritz values.')
               call svout (logfil, nev, bounds(np+1), ndigit, &
                    '_saup2: Ritz estimates of the "wanted" values ')
            end if
         end if

!
         if (ishift == 0) then
!
!           %-----------------------------------------------------%
!           | User specified shifts: reverse communication to     |
!           | compute the shifts. They are returned in the first  |
!           | NP locations of WORKL.                              |
!           %-----------------------------------------------------%
!
            ushift = .true.
            ido = 3
            go to 9000
         end if
!
   50    continue
!
!        %------------------------------------%
!        | Back from reverse communication;   |
!        | User specified shifts are returned |
!        | in WORKL(1:*NP)                   |
!        %------------------------------------%
!
         ushift = .false.
!
!
!        %---------------------------------------------------------%
!        | Move the NP shifts to the first NP locations of RITZ to |
!        | free up WORKL.  This is for the non-exact shift case;   |
!        | in the exact shift case, ssgets already handles this.   |
!        %---------------------------------------------------------%
!
         if (ishift == 0) call scopy (np, workl, 1, ritz, 1)
!
         if (msglvl > 2) then
            call ivout (logfil, 1, np, ndigit, &
                        '_saup2: The number of shifts to apply ')
            call svout (logfil, np, workl, ndigit, &
                        '_saup2: shifts selected')
            if (ishift == 1) then
               call svout (logfil, np, bounds, ndigit, &
                        '_saup2: corresponding Ritz estimates')
             end if
         end if
!
!        %---------------------------------------------------------%
!        | Apply the NP0 implicit shifts by QR bulge chasing.      |
!        | Each shift is applied to the entire tridiagonal matrix. |
!        | The first 2*N locations of WORKD are used as workspace. |
!        | After ssapps is done, we have a Lanczos                 |
!        | factorization of length NEV.                            |
!        %---------------------------------------------------------%
!
         call ssapps (n, nev, np, ritz, v, ldv, h, ldh, resid, q, ldq, &
              workd)
!
!        %---------------------------------------------%
!        | Compute the B-norm of the updated residual. |
!        | Keep B*RESID in WORKD(1:N) to be used in    |
!        | the first step of the next call to ssaitr.  |
!        %---------------------------------------------%
!
         cnorm = .true.
         call cpu_time (t2)
         if (bmat == 'G') then
            nbx = nbx + 1
            call scopy (n, resid, 1, workd(n+1), 1)
            ipntr(1) = n + 1
            ipntr(2) = 1
            ido = 2
!
!           %----------------------------------%
!           | Exit in order to compute B*RESID |
!           %----------------------------------%
!
            go to 9000
         else if (bmat == 'I') then
            call scopy (n, resid, 1, workd, 1)
         end if
!
  100    continue
!
!        %----------------------------------%
!        | Back from reverse communication; |
!        | WORKD(1:N) := B*RESID            |
!        %----------------------------------%
!
         if (bmat == 'G') then
            call cpu_time (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
!
         if (bmat == 'G') then
            rnorm = sdot (n, resid, 1, workd, 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat == 'I') then
            rnorm = snrm2(n, resid, 1)
         end if
         cnorm = .false.
  130    continue
!
         if (msglvl > 2) then
            call svout (logfil, 1, rnorm, ndigit, &
            '_saup2: B-norm of residual for NEV factorization')
            call svout (logfil, nev, h(1,2), ndigit, &
                 '_saup2: main diagonal of compressed H matrix')
            call svout (logfil, nev-1, h(2,1), ndigit, &
                 '_saup2: subdiagonal of compressed H matrix')
         end if
!
      go to 1000
!
!     %---------------------------------------------------------------%
!     |                                                               |
!     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
!     |                                                               |
!     %---------------------------------------------------------------%
!
 1100 continue
!
      mxiter = iter
      nev = nconv
!
 1200 continue
      ido = 99
!
!     %------------%
!     | Error exit |
!     %------------%
!
      call cpu_time (t1)
      tsaup2 = t1 - t0
!
 9000 continue
      return
!
!     %---------------%
!     | End of ssaup2 |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssaupd
!
!\Description:
!
!  Reverse communication interface for the Implicitly Restarted Arnoldi
!  Iteration.  For symmetric problems this reduces to a variant of the Lanczos
!  method.  This method has been designed to compute approximations to a
!  few eigenpairs of a linear operator OP that is real and symmetric
!  with respect to a real positive semi-definite symmetric matrix B,
!  i.e.
!
!       B*OP = (OP`)*B.
!
!  Another way to express this condition is
!
!       < x,OPy > = < OPx,y >  where < z,w > = z`Bw  .
!
!  In the standard eigenproblem B is the identity matrix.
!  ( A` denotes transpose of A)
!
!  The computed approximate eigenvalues are called Ritz values and
!  the corresponding approximate eigenvectors are called Ritz vectors.
!
!  ssaupd is usually called iteratively to solve one of the
!  following problems:
!
!  Mode 1:  A*x = lambda*x, A symmetric
!           ===> OP = A  and  B = I.
!
!  Mode 2:  A*x = lambda*M*x, A symmetric, M symmetric positive definite
!           ===> OP = inv[M]*A  and  B = M.
!           ===> (If M can be factored see remark 3 below)
!
!  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
!           ===> OP = (inv[K - sigma*M])*M  and  B = M.
!           ===> Shift-and-Invert mode
!
!  Mode 4:  K*x = lambda*KG*x, K symmetric positive semi-definite,
!           KG symmetric indefinite
!           ===> OP = (inv[K - sigma*KG])*K  and  B = K.
!           ===> Buckling mode
!
!  Mode 5:  A*x = lambda*M*x, A symmetric, M symmetric positive semi-definite
!           ===> OP = inv[A - sigma*M]*[A + sigma*M]  and  B = M.
!           ===> Cayley transformed mode
!
!  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
!        should be accomplished either by a direct method
!        using a sparse matrix factorization and solving
!
!           [A - sigma*M]*w = v  or M*w = v,
!
!        or through an iterative method for solving these
!        systems.  If an iterative method is used, the
!        convergence test must be more stringent than
!        the accuracy requirements for the eigenvalue
!        approximations.
!
!\Usage:
!  call ssaupd
!     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
!       IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!\Arguments
!  IDO     Integer.  (INPUT/OUTPUT)
!          Reverse communication flag.  IDO must be zero on the first
!          call to ssaupd.  IDO will be set internally to
!          indicate the type of operation to be performed.  Control is
!          then given back to the calling routine which has the
!          responsibility to carry out the requested operation and call
!          ssaupd with the result.  The operand is given in
!          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
!          (If Mode = 2 see remark 5 below)
!          -------------------------------------------------------------
!          IDO =  0: first call to the reverse communication interface
!          IDO = -1: compute  Y = OP * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    This is for the initialization phase to force the
!                    starting vector into the range of OP.
!          IDO =  1: compute  Y = OP * X where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!                    In mode 3,4 and 5, the vector B * X is already
!                    available in WORKD(ipntr(3)).  It does not
!                    need to be recomputed in forming OP * X.
!          IDO =  2: compute  Y = B * X  where
!                    IPNTR(1) is the pointer into WORKD for X,
!                    IPNTR(2) is the pointer into WORKD for Y.
!          IDO =  3: compute the IPARAM(8) shifts where
!                    IPNTR(11) is the pointer into WORKL for
!                    placing the shifts. See remark 6 below.
!          IDO = 99: done
!          -------------------------------------------------------------
!
!  BMAT    Character*1.  (INPUT)
!          BMAT specifies the type of the matrix B that defines the
!          semi-inner product for the operator OP.
!          B = 'I' -> standard eigenvalue problem A*x = lambda*x
!          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x
!
!  N       Integer.  (INPUT)
!          Dimension of the eigenproblem.
!
!  WHICH   Character*2.  (INPUT)
!          Specify which of the Ritz values of OP to compute.
!
!          'LA' - compute the NEV largest (algebraic) eigenvalues.
!          'SA' - compute the NEV smallest (algebraic) eigenvalues.
!          'LM' - compute the NEV largest (in magnitude) eigenvalues.
!          'SM' - compute the NEV smallest (in magnitude) eigenvalues.
!          'BE' - compute NEV eigenvalues, half from each end of the
!                 spectrum.  When NEV is odd, compute one more from the
!                 high end than from the low end.
!           (see remark 1 below)
!
!  NEV     Integer.  (INPUT)
!          Number of eigenvalues of OP to be computed. 0 < NEV < N.
!
!  TOL     Real  scalar.  (INPUT)
!          Stopping criterion: the relative accuracy of the Ritz value
!          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I)).
!          If TOL .LE. 0. is passed a default is set:
!          DEFAULT = SLAMCH('EPS')  (machine precision as computed
!                    by the LAPACK auxiliary subroutine SLAMCH).
!
!  RESID   Real  array of length N.  (INPUT/OUTPUT)
!          On INPUT:
!          If INFO .EQ. 0, a random initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          On OUTPUT:
!          RESID contains the final residual vector.
!
!  NCV     Integer.  (INPUT)
!          Number of columns of the matrix V (less than or equal to N).
!          This will indicate how many Lanczos vectors are generated
!          at each iteration.  After the startup phase in which NEV
!          Lanczos vectors are generated, the algorithm generates
!          NCV-NEV Lanczos vectors at each subsequent update iteration.
!          Most of the cost in generating each Lanczos vector is in the
!          matrix-vector product OP*x. (See remark 4 below).
!
!  V       Real  N by NCV array.  (OUTPUT)
!          The NCV columns of V contain the Lanczos basis vectors.
!
!  LDV     Integer.  (INPUT)
!          Leading dimension of V exactly as declared in the calling
!          program.
!
!  IPARAM  Integer array of length 11.  (INPUT/OUTPUT)
!          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
!          The shifts selected at each iteration are used to restart
!          the Arnoldi iteration in an implicit fashion.
!          -------------------------------------------------------------
!          ISHIFT = 0: the shifts are provided by the user via
!                      reverse communication.  The NCV eigenvalues of
!                      the current tridiagonal matrix T are returned in
!                      the part of WORKL array corresponding to RITZ.
!                      See remark 6 below.
!          ISHIFT = 1: exact shifts with respect to the reduced
!                      tridiagonal matrix T.  This is equivalent to
!                      restarting the iteration with a starting vector
!                      that is a linear combination of Ritz vectors
!                      associated with the "wanted" Ritz values.
!          -------------------------------------------------------------
!
!          IPARAM(2) = LEVEC
!          No longer referenced. See remark 2 below.
!
!          IPARAM(3) = MXITER
!          On INPUT:  maximum number of Arnoldi update iterations allowed.
!          On OUTPUT: actual number of Arnoldi update iterations taken.
!
!          IPARAM(4) = NB: blocksize to be used in the recurrence.
!          The code currently works only for NB = 1.
!
!          IPARAM(5) = NCONV: number of "converged" Ritz values.
!          This represents the number of Ritz values that satisfy
!          the convergence criterion.
!
!          IPARAM(6) = IUPD
!          No longer referenced. Implicit restarting is ALWAYS used.
!
!          IPARAM(7) = MODE
!          On INPUT determines what type of eigenproblem is being solved.
!          Must be 1,2,3,4,5; See under \Description of ssaupd for the
!          five modes available.
!
!          IPARAM(8) = NP
!          When ido = 3 and the user provides shifts through reverse
!          communication (IPARAM(1)=0), ssaupd returns NP, the number
!          of shifts the user is to provide. 0 < NP <=NCV-NEV. See Remark
!          6 below.
!
!          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
!          OUTPUT: NUMOP  = total number of OP*x operations,
!                  NUMOPB = total number of B*x operations if BMAT='G',
!                  NUMREO = total number of steps of re-orthogonalization.
!
!  IPNTR   Integer array of length 11.  (OUTPUT)
!          Pointer to mark the starting locations in the WORKD and WORKL
!          arrays for matrices/vectors used by the Lanczos iteration.
!          -------------------------------------------------------------
!          IPNTR(1): pointer to the current operand vector X in WORKD.
!          IPNTR(2): pointer to the current result vector Y in WORKD.
!          IPNTR(3): pointer to the vector B * X in WORKD when used in
!                    the shift-and-invert mode.
!          IPNTR(4): pointer to the next available location in WORKL
!                    that is untouched by the program.
!          IPNTR(5): pointer to the NCV by 2 tridiagonal matrix T in WORKL.
!          IPNTR(6): pointer to the NCV RITZ values array in WORKL.
!          IPNTR(7): pointer to the Ritz estimates in array WORKL associated
!                    with the Ritz values located in RITZ in WORKL.
!          IPNTR(11): pointer to the NP shifts in WORKL. See Remark 6 below.
!
!          Note: IPNTR(8:10) is only referenced by sseupd. See Remark 2.
!          IPNTR(8): pointer to the NCV RITZ values of the original system.
!          IPNTR(9): pointer to the NCV corresponding error bounds.
!          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
!                     of the tridiagonal matrix T. Only referenced by
!                     sseupd if RVEC = .TRUE. See Remarks.
!          -------------------------------------------------------------
!
!  WORKD   Real  work array of length 3*N.  (REVERSE COMMUNICATION)
!          Distributed array to be used in the basic Arnoldi iteration
!          for reverse communication.  The user should not use WORKD
!          as temporary workspace during the iteration. Upon termination
!          WORKD(1:N) contains B*RESID(1:N). If the Ritz vectors are desired
!          subroutine sseupd uses this output.
!          See Data Distribution Note below.
!
!  WORKL   Real  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.  See Data Distribution Note below.
!
!  LWORKL  Integer.  (INPUT)
!          LWORKL must be at least NCV**2 + 8*NCV .
!
!  INFO    Integer.  (INPUT/OUTPUT)
!          If INFO .EQ. 0, a randomly initial residual vector is used.
!          If INFO .NE. 0, RESID contains the initial residual vector,
!                          possibly from a previous run.
!          Error flag on output.
!          =  0: Normal exit.
!          =  1: Maximum number of iterations taken.
!                All possible eigenvalues of OP has been found. IPARAM(5)
!                returns the number of wanted converged Ritz values.
!          =  2: No longer an informational error. Deprecated starting
!                with release 2 of ARPACK.
!          =  3: No shifts could be applied during a cycle of the
!                Implicitly restarted Arnoldi iteration. One possibility
!                is to increase the size of NCV relative to NEV.
!                See remark 4 below.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV must be greater than NEV and less than or equal to N.
!          = -4: The maximum number of Arnoldi update iterations allowed
!                must be greater than zero.
!          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work array WORKL is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Informatinal error from LAPACK routine ssteqr.
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
!          = -12: IPARAM(1) must be equal to 0 or 1.
!          = -13: NEV and WHICH = 'BE' are incompatable.
!          = -9999: Could not build an Arnoldi factorization.
!                   IPARAM(5) returns the size of the current Arnoldi
!                   factorization. The user is advised to check that
!                   enough workspace and array storage has been allocated.
!
!
!\Remarks
!  1. The converged Ritz values are always returned in ascending
!     algebraic order.  The computed Ritz values are approximate
!     eigenvalues of OP.  The selection of WHICH should be made
!     with this in mind when Mode = 3,4,5.  After convergence,
!     approximate eigenvalues of the original problem may be obtained
!     with the ARPACK subroutine sseupd.
!
!  2. If the Ritz vectors corresponding to the converged Ritz values
!     are needed, the user must call sseupd immediately following completion
!     of ssaupd. This is new starting with version 2.1 of ARPACK.
!
!  3. If M can be factored into a Cholesky factorization M = LL`
!     then Mode = 2 should not be selected.  Instead one should use
!     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
!     linear systems should be solved with L and L` rather
!     than computing inverses.  After convergence, an approximate
!     eigenvector z of the original problem is recovered by solving
!     L`z = x  where x is a Ritz vector of OP.
!
!  4. At present there is no a-priori analysis to guide the selection
!     of NCV relative to NEV.  The only formal requrement is that NCV > NEV.
!     However, it is recommended that NCV >= 2*NEV.  If many problems of
!     the same type are to be solved, one should experiment with increasing
!     NCV while keeping NEV fixed for a given test problem.  This will
!     usually decrease the required number of OP*x operations but it
!     also increases the work and storage required to maintain the orthogonal
!     basis vectors.   The optimal "cross-over" with respect to CPU time
!     is problem dependent and must be determined empirically.
!
!  5. If IPARAM(7) = 2 then in the Reverse commuication interface the user
!     must do the following. When IDO = 1, Y = OP * X is to be computed.
!     When IPARAM(7) = 2 OP = inv(B)*A. After computing A*X the user
!     must overwrite X with A*X. Y is then the solution to the linear set
!     of equations B*Y = A*X.
!
!  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
!     NP = IPARAM(8) shifts in locations:
!     1   WORKL(IPNTR(11))
!     2   WORKL(IPNTR(11)+1)
!                        .
!                        .
!                        .
!     NP  WORKL(IPNTR(11)+NP-1).
!
!     The eigenvalues of the current tridiagonal matrix are located in
!     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
!     order defined by WHICH. The associated Ritz estimates are located in
!     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
!
!-----------------------------------------------------------------------
!
!\Data Distribution Note:
!
!  Fortran-D syntax:
!  ================
!  REAL       RESID(N), V(LDV,NCV), WORKD(3*N), WORKL(LWORKL)
!  DECOMPOSE  D1(N), D2(N,NCV)
!  ALIGN      RESID(I) with D1(I)
!  ALIGN      V(I,J)   with D2(I,J)
!  ALIGN      WORKD(I) with D1(I)     range (1:N)
!  ALIGN      WORKD(I) with D1(I-N)   range (N+1:2*N)
!  ALIGN      WORKD(I) with D1(I-2*N) range (2*N+1:3*N)
!  DISTRIBUTE D1(BLOCK), D2(BLOCK,:)
!  REPLICATED WORKL(LWORKL)
!
!  Cray MPP syntax:
!  ===============
!  REAL       RESID(N), V(LDV,NCV), WORKD(N,3), WORKL(LWORKL)
!  SHARED     RESID(BLOCK), V(BLOCK,:), WORKD(BLOCK,:)
!  REPLICATED WORKL(LWORKL)
!
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!  8. R.B. Lehoucq, D.C. Sorensen, "Implementation of Some Spectral
!     Transformations in a k-Step Arnoldi Method". In Preparation.
!
!\Routines called:
!     ssaup2  ARPACK routine that implements the Implicitly Restarted
!             Arnoldi Iteration.
!     sstats  ARPACK routine that initialize timing and other statistics
!             variables.
!     ivout   ARPACK utility routine that prints integers.
!     svout   ARPACK utility routine that prints vectors.
!     slamch  LAPACK routine that determines machine constants.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/15/93: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: saupd.F   SID: 2.8   DATE OF SID: 04/10/01   RELEASE: 2
!
!\Remarks
!     1. None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssaupd &
         ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
           ipntr, workd, workl, lworkl, info )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat*1, which*2
      integer    ido, info, ldv, lworkl, n, ncv, nev
      Real &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(11), ipntr(11)
      Real &
                 resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0 , zero = 0.0E+0 )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    bounds, ierr, ih, iq, ishift, iupd, iw, &
                 ldh, ldq, msglvl, mxiter, mode, nb, &
                 nev0, next, np, ritz, j
      save       bounds, ierr, ih, iq, ishift, iupd, iw, &
                 ldh, ldq, msglvl, mxiter, mode, nb, &
                 nev0, next, np, ritz
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   ssaup2,  svout, ivout, sstats
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slamch
      external   slamch
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      if (ido == 0) then
!
!        %-------------------------------%
!        | Initialize timing statistics  |
!        | & message level for debugging |
!        %-------------------------------%
!
         call sstats
         call cpu_time (t0)
         msglvl = msaupd
!
         ierr   = 0
         ishift = iparam(1)
         mxiter = iparam(3)
!         nb     = iparam(4)
         nb     = 1
!
!        %--------------------------------------------%
!        | Revision 2 performs only implicit restart. |
!        %--------------------------------------------%
!
         iupd   = 1
         mode   = iparam(7)
!
!        %----------------%
!        | Error checking |
!        %----------------%
!
         if (n <= 0) then
            ierr = -1
         else if (nev <= 0) then
            ierr = -2
         else if (ncv <= nev .or.  ncv > n) then
            ierr = -3
         end if
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        %----------------------------------------------%
!
         np     = ncv - nev
!
         if (mxiter <= 0)                     ierr = -4
         if (which /= 'LM' .and. &
             which /= 'SM' .and. &
             which /= 'LA' .and. &
             which /= 'SA' .and. &
             which /= 'BE')                   ierr = -5
         if (bmat /= 'I' .and. bmat .ne. 'G') ierr = -6
!
         if (lworkl < ncv**2 + 8*ncv)        ierr = -7
         if (mode < 1 .or. mode > 5) then
                                                ierr = -10
         else if (mode == 1 .and. bmat .eq. 'G') then
                                                ierr = -11
         else if (ishift < 0 .or. ishift > 1) then
                                                ierr = -12
         else if (nev == 1 .and. which .eq. 'BE') then
                                                ierr = -13
         end if
!
!        %------------%
!        | Error Exit |
!        %------------%
!
         if (ierr /= 0) then
            info = ierr
            ido  = 99
            go to 9000
         end if
!
!        %------------------------%
!        | Set default parameters |
!        %------------------------%
!
         if (nb <= 0)                         nb = 1
         if (tol <= zero)                     tol = slamch('EpsMach')
!
!        %----------------------------------------------%
!        | NP is the number of additional steps to      |
!        | extend the length NEV Lanczos factorization. |
!        | NEV0 is the local variable designating the   |
!        | size of the invariant subspace desired.      |
!        %----------------------------------------------%
!
         np     = ncv - nev
         nev0   = nev
!
!        %-----------------------------%
!        | Zero out internal workspace |
!        %-----------------------------%
!
         do j = 1, ncv**2 + 8*ncv
            workl(j) = zero
         end do
!
!        %-------------------------------------------------------%
!        | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
!        | etc... and the remaining workspace.                   |
!        | Also update pointer to be used on output.             |
!        | Memory is laid out as follows:                        |
!        | workl(1:2*ncv) := generated tridiagonal matrix        |
!        | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
!        | workl(3*ncv+1:3*ncv+ncv) := computed error bounds     |
!        | workl(4*ncv+1:4*ncv+ncv*ncv) := rotation matrix Q     |
!        | workl(4*ncv+ncv*ncv+1:7*ncv+ncv*ncv) := workspace     |
!        %-------------------------------------------------------%
!
         ldh    = ncv
         ldq    = ncv
         ih     = 1
         ritz   = ih     + 2*ldh
         bounds = ritz   + ncv
         iq     = bounds + ncv
         iw     = iq     + ncv**2
         next   = iw     + 3*ncv
!
         ipntr(4) = next
         ipntr(5) = ih
         ipntr(6) = ritz
         ipntr(7) = bounds
         ipntr(11) = iw
      end if
!
!     %-------------------------------------------------------%
!     | Carry out the Implicitly restarted Lanczos Iteration. |
!     %-------------------------------------------------------%
!
      call ssaup2 &
         ( ido, bmat, n, which, nev0, np, tol, resid, mode, iupd, &
           ishift, mxiter, v, ldv, workl(ih), ldh, workl(ritz), &
           workl(bounds), workl(iq), ldq, workl(iw), ipntr, workd, &
           info )
!
!     %--------------------------------------------------%
!     | ido /= 99 implies use of reverse communication |
!     | to compute operations involving OP or shifts.    |
!     %--------------------------------------------------%
!
      if (ido == 3) iparam(8) = np
      if (ido /= 99) go to 9000
!
      iparam(3) = mxiter
      iparam(5) = np
      iparam(9) = nopx
      iparam(10) = nbx
      iparam(11) = nrorth
!
!     %------------------------------------%
!     | Exit if there was an informational |
!     | error within ssaup2.               |
!     %------------------------------------%
!
      if (info < 0) go to 9000
      if (info == 2) info = 3
!
      if (msglvl > 0) then
         call ivout (logfil, 1, mxiter, ndigit, &
                     '_saupd: number of update iterations taken')
         call ivout (logfil, 1, np, ndigit, &
                     '_saupd: number of "converged" Ritz values')
         call svout (logfil, np, workl(Ritz), ndigit, &
                     '_saupd: final Ritz values')
         call svout (logfil, np, workl(Bounds), ndigit, &
                     '_saupd: corresponding error bounds')
      end if
!
      call cpu_time (t1)
      tsaupd = t1 - t0
!
      if (msglvl > 0) then
!
!        %--------------------------------------------------------%
!        | Version Number & Version Date are defined in version.h |
!        %--------------------------------------------------------%
!
         write (6,1000)
         write (6,1100) mxiter, nopx, nbx, nrorth, nitref, nrstrt, &
                        tmvopx, tmvbx, tsaupd, tsaup2, tsaitr, titref, &
                        tgetv0, tseigt, tsgets, tsapps, tsconv
 1000    format (//, &
            5x, '==========================================',/ &
            5x, '= Symmetric implicit Arnoldi update code =',/ &
            5x, '= Version Number:', ' 2.4' , 19x, ' =',/ &
            5x, '= Version Date:  ', ' 07/31/96' , 14x, ' =',/ &
            5x, '==========================================',/ &
            5x, '= Summary of timing statistics           =',/ &
            5x, '==========================================',//)
 1100    format ( &
            5x, 'Total number update iterations             = ', i5,/ &
            5x, 'Total number of OP*x operations            = ', i5,/ &
            5x, 'Total number of B*x operations             = ', i5,/ &
            5x, 'Total number of reorthogonalization steps  = ', i5,/ &
            5x, 'Total number of iterative refinement steps = ', i5,/ &
            5x, 'Total number of restart steps              = ', i5,/ &
            5x, 'Total time in user OP*x operation          = ', f12.6,/ &
            5x, 'Total time in user B*x operation           = ', f12.6,/ &
            5x, 'Total time in Arnoldi update routine       = ', f12.6,/ &
            5x, 'Total time in saup2 routine                = ', f12.6,/ &
            5x, 'Total time in basic Arnoldi iteration loop = ', f12.6,/ &
            5x, 'Total time in reorthogonalization phase    = ', f12.6,/ &
            5x, 'Total time in (re)start vector generation  = ', f12.6,/ &
            5x, 'Total time in trid eigenvalue subproblem   = ', f12.6,/ &
            5x, 'Total time in getting the shifts           = ', f12.6,/ &
            5x, 'Total time in applying the shifts          = ', f12.6,/ &
            5x, 'Total time in convergence testing          = ', f12.6)
      end if
!
 9000 continue
!
      return
!
!     %---------------%
!     | End of ssaupd |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssconv
!
!\Description:
!  Convergence testing for the symmetric Arnoldi eigenvalue routine.
!
!\Usage:
!  call ssconv
!     ( N, RITZ, BOUNDS, TOL, NCONV )
!
!\Arguments
!  N       Integer.  (INPUT)
!          Number of Ritz values to check for convergence.
!
!  RITZ    Real array of length N.  (INPUT)
!          The Ritz values to be checked for convergence.
!
!  BOUNDS  Real array of length N.  (INPUT)
!          Ritz estimates associated with the Ritz values in RITZ.
!
!  TOL     Real scalar.  (INPUT)
!          Desired relative accuracy for a Ritz value to be considered
!          "converged".
!
!  NCONV   Integer scalar.  (OUTPUT)
!          Number of "converged" Ritz values.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines called:
!     slamch  LAPACK routine that determines machine constants.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.4, this routine no longer uses the
!        Parlett strategy using the gap conditions.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssconv (n, ritz, bounds, tol, nconv)
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    n, nconv
      Real &
                 tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 ritz(n), bounds(n)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i
      Real &
                 temp, eps23
!
!     %-------------------%
!     | External routines |
!     %-------------------%
!
      Real &
                 slamch
      external   slamch

!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    abs
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      call cpu_time (t0)

      eps23 = slamch('Epsilon-Machine')
      eps23 = eps23**(2.0E+0 / 3.0E+0)

      nconv  = 0
      do i = 1, n
!
!        %-----------------------------------------------------%
!        | The i-th Ritz value is considered "converged"       |
!        | when: bounds(i) <= TOL*max(eps23, abs(ritz(i)))   |
!        %-----------------------------------------------------%
!
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) <= tol*temp ) then
            nconv = nconv + 1
         end if

      end do
!
      call cpu_time (t1)
      tsconv = tsconv + (t1 - t0)
!
      return
!
!     %---------------%
!     | End of ssconv |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: sseigt
!
!\Description:
!  Compute the eigenvalues of the current symmetric tridiagonal matrix
!  and the corresponding error bounds given the current residual norm.
!
!\Usage:
!  call sseigt
!     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
!
!\Arguments
!  RNORM   Real scalar.  (INPUT)
!          RNORM contains the residual norm corresponding to the current
!          symmetric tridiagonal matrix H.
!
!  N       Integer.  (INPUT)
!          Size of the symmetric tridiagonal matrix H.
!
!  H       Real N by 2 array.  (INPUT)
!          H contains the symmetric tridiagonal matrix with the
!          subdiagonal in the first column starting at H(2,1) and the
!          main diagonal in second column.
!
!  LDH     Integer.  (INPUT)
!          Leading dimension of H exactly as declared in the calling
!          program.
!
!  EIG     Real array of length N.  (OUTPUT)
!          On output, EIG contains the N eigenvalues of H possibly
!          unsorted.  The BOUNDS arrays are returned in the
!          same sorted order as EIG.
!
!  BOUNDS  Real array of length N.  (OUTPUT)
!          On output, BOUNDS contains the error estimates corresponding
!          to the eigenvalues EIG.  This is equal to RNORM times the
!          last components of the eigenvectors corresponding to the
!          eigenvalues in EIG.
!
!  WORKL   Real work array of length 3*N.  (WORKSPACE)
!          Private (replicated) array on each PE or array allocated on
!          the front end.
!
!  IERR    Integer.  (OUTPUT)
!          Error exit flag from sstqrb.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     sstqrb  ARPACK routine that computes the eigenvalues and the
!             last components of the eigenvectors of a symmetric
!             and tridiagonal matrix.
!     svout   ARPACK utility routine that prints vectors.
!     scopy   Level 1 BLAS that copies one vector to another.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.4'
!
!\SCCS Information: @(#)
! FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     None
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine sseigt &
         ( rnorm, n, h, ldh, eig, bounds, workl, ierr )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    ierr, ldh, n
      Real &
                 rnorm
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 eig(n), bounds(n), h(ldh,2), workl(3*n)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 zero
      parameter (zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    k, msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy, sstqrb, svout
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = mseigt
!
      if (msglvl > 0) then
         call svout (logfil, n, h(1,2), ndigit, &
                    '_seigt: main diagonal of matrix H')
         if (n > 1) then
         call svout (logfil, n-1, h(2,1), ndigit, &
                    '_seigt: sub diagonal of matrix H')
         end if
      end if
!
      call scopy  (n, h(1,2), 1, eig, 1)
      call scopy  (n-1, h(2,1), 1, workl, 1)
      call sstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr /= 0) go to 9000
      if (msglvl > 1) then
         call svout (logfil, n, bounds, ndigit, &
                    '_seigt: last row of the eigenvector matrix for H')
      end if
!
!     %-----------------------------------------------%
!     | Finally determine the error bounds associated |
!     | with the n Ritz values of H.                  |
!     %-----------------------------------------------%
!
      do 30 k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
   30 continue
!
      call cpu_time (t1)
      tseigt = tseigt + (t1 - t0)
!
 9000 continue
      return
!
!     %---------------%
!     | End of sseigt |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssesrt
!
!\Description:
!  Sort the array X in the order specified by WHICH and optionally
!  apply the permutation to the columns of the matrix A.
!
!\Usage:
!  call ssesrt
!     ( WHICH, APPLY, N, X, NA, A, LDA)
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X is sorted into increasing order of magnitude.
!          'SM' -> X is sorted into decreasing order of magnitude.
!          'LA' -> X is sorted into increasing order of algebraic.
!          'SA' -> X is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to A.
!          APPLY = .FALSE. -> do not apply the sorted order to A.
!
!  N       Integer.  (INPUT)
!          Dimension of the array X.
!
!  X      Real array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  NA      Integer.  (INPUT)
!          Number of rows of the matrix A.
!
!  A      Real array of length NA by N.  (INPUT/OUTPUT)
!
!  LDA     Integer.  (INPUT)
!          Leading dimension of A.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Routines
!     sswap  Level 1 BLAS that swaps the contents of two vectors.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/15/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO and
!               the ARPACK code ssortr
!
!\SCCS Information: @(#)
! FILE: sesrt.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssesrt (which, apply, n, x, na, a, lda)
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      logical    apply
      integer    lda, n, na
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 x(0:n-1), a(lda, 0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, igap, j
      Real &
                 temp
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   sswap
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      igap = n / 2
!
      if (which == 'SA') then
!
!        X is sorted into decreasing order of algebraic.
!
   10    continue
         if (igap == 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if (j<0) go to 30
!
            if (x(j)<x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
!
      else if (which == 'SM') then
!
!        X is sorted into decreasing order of magnitude.
!
   40    continue
         if (igap == 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j<0) go to 60
!
            if (abs(x(j))<abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
!
      else if (which == 'LA') then
!
!        X is sorted into increasing order of algebraic.
!
   70    continue
         if (igap == 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j<0) go to 90
!
            if (x(j)>x(j+igap)) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
!
      else if (which == 'LM') then
!
!        X is sorted into increasing order of magnitude.
!
  100    continue
         if (igap == 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j<0) go to 120
!
            if (abs(x(j))>abs(x(j+igap))) then
               temp = x(j)
               x(j) = x(j+igap)
               x(j+igap) = temp
               if (apply) call sswap( na, a(1, j), 1, a(1,j+igap), 1)
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
!
 9000 continue
      return
!
!     %---------------%
!     | End of ssesrt |
!     %---------------%
!
      end
!\BeginDoc
!
!\Name: sseupd
!
!\Description:
!
!  This subroutine returns the converged approximations to eigenvalues
!  of A*z = lambda*B*z and (optionally):
!
!      (1) the corresponding approximate eigenvectors,
!
!      (2) an orthonormal (Lanczos) basis for the associated approximate
!          invariant subspace,
!
!      (3) Both.
!
!  There is negligible additional cost to obtain eigenvectors.  An orthonormal
!  (Lanczos) basis is always computed.  There is an additional storage cost
!  of n*nev if both are requested (in this case a separate array Z must be
!  supplied).
!
!  These quantities are obtained from the Lanczos factorization computed
!  by SSAUPD for the linear operator OP prescribed by the MODE selection
!  (see IPARAM(7) in SSAUPD documentation.)  SSAUPD must be called before
!  this routine is called. These approximate eigenvalues and vectors are
!  commonly called Ritz values and Ritz vectors respectively.  They are
!  referred to as such in the comments that follow.   The computed orthonormal
!  basis for the invariant subspace corresponding to these Ritz values is
!  referred to as a Lanczos basis.
!
!  See documentation in the header of the subroutine SSAUPD for a definition
!  of OP as well as other terms and the relation of computed Ritz values
!  and vectors of OP with respect to the given problem  A*z = lambda*B*z.
!
!  The approximate eigenvalues of the original problem are returned in
!  ascending algebraic order.  The user may elect to call this routine
!  once for each desired Ritz vector and store it peripherally if desired.
!  There is also the option of computing a selected set of these vectors
!  with a single call.
!
!\Usage:
!  call sseupd
!     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, BMAT, N, WHICH, NEV, TOL,
!       RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD, WORKL, LWORKL, INFO )
!
!  RVEC    LOGICAL  (INPUT)
!          Specifies whether Ritz vectors corresponding to the Ritz value
!          approximations to the eigenproblem A*z = lambda*B*z are computed.
!
!             RVEC = .FALSE.     Compute Ritz values only.
!
!             RVEC = .TRUE.      Compute Ritz vectors.
!
!  HOWMNY  Character*1  (INPUT)
!          Specifies how many Ritz vectors are wanted and the form of Z
!          the matrix of Ritz vectors. See remark 1 below.
!          = 'A': compute NEV Ritz vectors;
!          = 'S': compute some of the Ritz vectors, specified
!                 by the logical array SELECT.
!
!  SELECT  Logical array of dimension NCV.  (INPUT/WORKSPACE)
!          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
!          computed. To select the Ritz vector corresponding to a
!          Ritz value D(j), SELECT(j) must be set to .TRUE..
!          If HOWMNY = 'A' , SELECT is used as a workspace for
!          reordering the Ritz values.
!
!  D       Real  array of dimension NEV.  (OUTPUT)
!          On exit, D contains the Ritz value approximations to the
!          eigenvalues of A*z = lambda*B*z. The values are returned
!          in ascending order. If IPARAM(7) = 3,4,5 then D represents
!          the Ritz values of OP computed by ssaupd transformed to
!          those of the original eigensystem A*z = lambda*B*z. If
!          IPARAM(7) = 1,2 then the Ritz values of OP are the same
!          as the those of A*z = lambda*B*z.
!
!  Z       Real  N by NEV array if HOWMNY = 'A'.  (OUTPUT)
!          On exit, Z contains the B-orthonormal Ritz vectors of the
!          eigensystem A*z = lambda*B*z corresponding to the Ritz
!          value approximations.
!          If  RVEC = .FALSE. then Z is not referenced.
!          NOTE: The array Z may be set equal to first NEV columns of the
!          Arnoldi/Lanczos basis array V computed by SSAUPD.
!
!  LDZ     Integer.  (INPUT)
!          The leading dimension of the array Z.  If Ritz vectors are
!          desired, then  LDZ >=  max( 1, N ).  In any case,  LDZ .ge. 1.
!
!  SIGMA   Real   (INPUT)
!          If IPARAM(7) = 3,4,5 represents the shift. Not referenced if
!          IPARAM(7) = 1 or 2.
!
!
!  **** The remaining arguments MUST be the same as for the   ****
!  **** call to SSAUPD that was just completed.               ****
!
!  NOTE: The remaining arguments
!
!           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
!           WORKD, WORKL, LWORKL, INFO
!
!         must be passed directly to SSEUPD following the last call
!         to SSAUPD.  These arguments MUST NOT BE MODIFIED between
!         the the last call to SSAUPD and the call to SSEUPD.
!
!  Two of these parameters (WORKL, INFO) are also output parameters:
!
!  WORKL   Real  work array of length LWORKL.  (OUTPUT/WORKSPACE)
!          WORKL(1:4*ncv) contains information obtained in
!          ssaupd.  They are not changed by sseupd.
!          WORKL(4*ncv+1:ncv*ncv+8*ncv) holds the
!          untransformed Ritz values, the computed error estimates,
!          and the associated eigenvector matrix of H.
!
!          Note: IPNTR(8:10) contains the pointer into WORKL for addresses
!          of the above information computed by sseupd.
!          -------------------------------------------------------------
!          IPNTR(8): pointer to the NCV RITZ values of the original system.
!          IPNTR(9): pointer to the NCV corresponding error bounds.
!          IPNTR(10): pointer to the NCV by NCV matrix of eigenvectors
!                     of the tridiagonal matrix T. Only referenced by
!                     sseupd if RVEC = .TRUE. See Remarks.
!          -------------------------------------------------------------
!
!  INFO    Integer.  (OUTPUT)
!          Error flag on output.
!          =  0: Normal exit.
!          = -1: N must be positive.
!          = -2: NEV must be positive.
!          = -3: NCV must be greater than NEV and less than or equal to N.
!          = -5: WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.
!          = -6: BMAT must be one of 'I' or 'G'.
!          = -7: Length of private work WORKL array is not sufficient.
!          = -8: Error return from trid. eigenvalue calculation;
!                Information error from LAPACK routine ssteqr.
!          = -9: Starting vector is zero.
!          = -10: IPARAM(7) must be 1,2,3,4,5.
!          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
!          = -12: NEV and WHICH = 'BE' are incompatible.
!          = -14: SSAUPD did not find any eigenvalues to sufficient
!                 accuracy.
!          = -15: HOWMNY must be one of 'A' or 'S' if RVEC = .true.
!          = -16: HOWMNY = 'S' not yet implemented
!          = -17: SSEUPD got a different count of the number of converged
!                 Ritz values than SSAUPD got.  This indicates the user
!                 probably made an error in passing data from SSAUPD to
!                 SSEUPD or that the data was modified before entering
!                 SSEUPD.
!
!\BeginLib
!
!\References:
!  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
!     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
!     pp 357-385.
!  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
!     Restarted Arnoldi Iteration", Rice University Technical Report
!     TR95-13, Department of Computational and Applied Mathematics.
!  3. B.N. Parlett, "The Symmetric Eigenvalue Problem". Prentice-Hall,
!     1980.
!  4. B.N. Parlett, B. Nour-Omid, "Towards a Black Box Lanczos Program",
!     Computer Physics Communications, 53 (1989), pp 169-179.
!  5. B. Nour-Omid, B.N. Parlett, T. Ericson, P.S. Jensen, "How to
!     Implement the Spectral Transformation", Math. Comp., 48 (1987),
!     pp 663-673.
!  6. R.G. Grimes, J.G. Lewis and H.D. Simon, "A Shifted Block Lanczos
!     Algorithm for Solving Sparse Symmetric Generalized Eigenproblems",
!     SIAM J. Matr. Anal. Apps.,  January (1993).
!  7. L. Reichel, W.B. Gragg, "Algorithm 686: FORTRAN Subroutines
!     for Updating the QR decomposition", ACM TOMS, December 1990,
!     Volume 16 Number 4, pp 369-377.
!
!\Remarks
!  1. The converged Ritz values are always returned in increasing
!     (algebraic) order.
!
!  2. Currently only HOWMNY = 'A' is implemented. It is included at this
!     stage for the user who wants to incorporate it.
!
!\Routines called:
!     ssesrt  ARPACK routine that sorts an array X, and applies the
!             corresponding permutation to a matrix A.
!     ssortr  ssortr  ARPACK sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     svout   ARPACK utility routine that prints vectors.
!     sgeqr2  LAPACK routine that computes the QR factorization of
!             a matrix.
!     slacpy  LAPACK matrix copy routine.
!     slamch  LAPACK routine that determines machine constants.
!     sorm2r  LAPACK routine that applies an orthogonal matrix in
!             factored form.
!     ssteqr  LAPACK routine that computes eigenvalues and eigenvectors
!             of a tridiagonal matrix.
!     sger    Level 2 BLAS rank one update to a matrix.
!     scopy   Level 1 BLAS that copies one vector to another .
!     snrm2   Level 1 BLAS that computes the norm of a vector.
!     sscal   Level 1 BLAS that scales a vector.
!     sswap   Level 1 BLAS that swaps the contents of two vectors.

!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Chao Yang                    Houston, Texas
!     Dept. of Computational &
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/15/93: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: seupd.F   SID: 2.11   DATE OF SID: 04/10/01   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
      subroutine sseupd(rvec  , howmny, select, d    , &
                         z     , ldz   , sigma , bmat , &
                         n     , which , nev   , tol  , &
                         resid , ncv   , v     , ldv  , &
                         iparam, ipntr , workd , workl, &
                         lworkl, info )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character  bmat, howmny, which*2
      logical    rvec
      integer    info, ldz, ldv, lworkl, n, ncv, nev
      Real &
                 sigma, tol
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      integer    iparam(7), ipntr(11)
      logical    select(ncv)
      Real &
                 d(nev)     , resid(n)  , v(ldv,ncv), &
                 z(ldz, nev), workd(2*n), workl(lworkl)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0 , zero = 0.0E+0 )
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      character  type*6
      integer    bounds , ierr   , ih    , ihb   , ihd   , &
                 iq     , iw     , j     , k     , ldh   , &
                 ldq    , mode   , msglvl, nconv , next  , &
                 ritz   , irz    , ibd   , np    , ishift, &
                 leftptr, rghtptr, numcnv, jj
      Real &
                 bnorm2 , rnorm, temp, temp1, eps23
      logical    reord
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   scopy , sger  , sgeqr2, slacpy, sorm2r, sscal, &
                 ssesrt, ssteqr, sswap , svout , ivout , ssortr
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 snrm2, slamch
      external   snrm2, slamch
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %------------------------%
!     | Set default parameters |
!     %------------------------%
!
      msglvl = mseupd
      mode = iparam(7)
      nconv = iparam(5)
      info = 0
!
!     %--------------%
!     | Quick return |
!     %--------------%
!
      if (nconv == 0) go to 9000
      ierr = 0
!
      if (nconv <= 0)                        ierr = -14
      if (n <= 0)                            ierr = -1
      if (nev <= 0)                          ierr = -2
      if (ncv <= nev .or.  ncv > n)       ierr = -3
      if (which /= 'LM' .and. &
          which /= 'SM' .and. &
          which /= 'LA' .and. &
          which /= 'SA' .and. &
          which /= 'BE')                     ierr = -5
      if (bmat /= 'I' .and. bmat .ne. 'G')   ierr = -6
      if ( (howmny /= 'A' .and. &
                 howmny /= 'P' .and. &
                 howmny /= 'S') .and. rvec ) &
                                               ierr = -15
      if (rvec .and. howmny == 'S')           ierr = -16
!
      if (rvec .and. lworkl < ncv**2+8*ncv) ierr = -7
!
      if (mode == 1 .or. mode .eq. 2) then
         type = 'REGULR'
      else if (mode == 3 ) then
         type = 'SHIFTI'
      else if (mode == 4 ) then
         type = 'BUCKLE'
      else if (mode == 5 ) then
         type = 'CAYLEY'
      else
                                               ierr = -10
      end if
      if (mode == 1 .and. bmat .eq. 'G')     ierr = -11
      if (nev == 1 .and. which .eq. 'BE')    ierr = -12
!
!     %------------%
!     | Error Exit |
!     %------------%
!
      if (ierr /= 0) then
         info = ierr
         go to 9000
      end if
!
!     %-------------------------------------------------------%
!     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q  |
!     | etc... and the remaining workspace.                   |
!     | Also update pointer to be used on output.             |
!     | Memory is laid out as follows:                        |
!     | workl(1:2*ncv) := generated tridiagonal matrix H      |
!     |       The subdiagonal is stored in workl(2:ncv).      |
!     |       The dead spot is workl(1) but upon exiting      |
!     |       ssaupd stores the B-norm of the last residual   |
!     |       vector in workl(1). We use this !!!             |
!     | workl(2*ncv+1:2*ncv+ncv) := ritz values               |
!     |       The wanted values are in the first NCONV spots. |
!     | workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates   |
!     |       The wanted values are in the first NCONV spots. |
!     | NOTE: workl(1:4*ncv) is set by ssaupd and is not      |
!     |       modified by sseupd.                             |
!     %-------------------------------------------------------%
!
!     %-------------------------------------------------------%
!     | The following is used and set by sseupd.              |
!     | workl(4*ncv+1:4*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the diagonal of H. Upon EXIT contains the NCV   |
!     |       Ritz values of the original system. The first   |
!     |       NCONV spots have the wanted values. If MODE =   |
!     |       1 or 2 then will equal workl(2*ncv+1:3*ncv).    |
!     | workl(5*ncv+1:5*ncv+ncv) := used as workspace during  |
!     |       computation of the eigenvectors of H. Stores    |
!     |       the subdiagonal of H. Upon EXIT contains the    |
!     |       NCV corresponding Ritz estimates of the         |
!     |       original system. The first NCONV spots have the |
!     |       wanted values. If MODE = 1,2 then will equal    |
!     |       workl(3*ncv+1:4*ncv).                           |
!     | workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is  |
!     |       the eigenvector matrix for H as returned by     |
!     |       ssteqr. Not referenced if RVEC = .False.        |
!     |       Ordering follows that of workl(4*ncv+1:5*ncv)   |
!     | workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=         |
!     |       Workspace. Needed by ssteqr and by sseupd.      |
!     | GRAND total of NCV*(NCV+8) locations.                 |
!     %-------------------------------------------------------%
!
!
      ih     = ipntr(5)
      ritz   = ipntr(6)
      bounds = ipntr(7)
      ldh    = ncv
      ldq    = ncv
      ihd    = bounds + ldh
      ihb    = ihd    + ldh
      iq     = ihb    + ldh
      iw     = iq     + ldh*ncv
      next   = iw     + 2*ncv
      ipntr(4)  = next
      ipntr(8)  = ihd
      ipntr(9)  = ihb
      ipntr(10) = iq
!
!     %----------------------------------------%
!     | irz points to the Ritz values computed |
!     |     by _seigt before exiting _saup2.   |
!     | ibd points to the Ritz estimates       |
!     |     computed by _seigt before exiting  |
!     |     _saup2.                            |
!     %----------------------------------------%
!
      irz = ipntr(11)+ncv
      ibd = irz+ncv
!
!
!     %---------------------------------%
!     | Set machine dependent constant. |
!     %---------------------------------%
!
      eps23 = slamch('Epsilon-Machine')
      eps23 = eps23**(2.0E+0  / 3.0E+0 )
!
!     %---------------------------------------%
!     | RNORM is B-norm of the RESID(1:N).    |
!     | BNORM2 is the 2 norm of B*RESID(1:N). |
!     | Upon exit of ssaupd WORKD(1:N) has    |
!     | B*RESID(1:N).                         |
!     %---------------------------------------%
!
      rnorm = workl(ih)
      if (bmat == 'I') then
         bnorm2 = rnorm
      else if (bmat == 'G') then
         bnorm2 = snrm2(n, workd, 1)
      end if

      if (msglvl > 2) then
         call svout(logfil, ncv, workl(irz), ndigit, &
         '_seupd: Ritz values passed in from _SAUPD.')
         call svout(logfil, ncv, workl(ibd), ndigit, &
         '_seupd: Ritz estimates passed in from _SAUPD.')
      end if

      if (rvec) then

         reord = .false.
!
!        %---------------------------------------------------%
!        | Use the temporary bounds array to store indices   |
!        | These will be used to mark the select array later |
!        %---------------------------------------------------%
!
         do j = 1,ncv
            workl(bounds+j-1) = j
            select(j) = .false.
         end do
!
!        %-------------------------------------%
!        | Select the wanted Ritz values.      |
!        | Sort the Ritz values so that the    |
!        | wanted ones appear at the tailing   |
!        | NEV positions of workl(irr) and     |
!        | workl(iri).  Move the corresponding |
!        | error estimates in workl(bound)     |
!        | accordingly.                        |
!        %-------------------------------------%
!
         np     = ncv - nev
         ishift = 0
         call ssgets(ishift, which       , nev          , &
                      np    , workl(irz)  , workl(bounds), &
                      workl)
!
         if (msglvl > 2) then
            call svout(logfil, ncv, workl(irz), ndigit, &
            '_seupd: Ritz values after calling _SGETS.')
            call svout(logfil, ncv, workl(bounds), ndigit, &
            '_seupd: Ritz value indices after calling _SGETS.')
         end if
!
!        %-----------------------------------------------------%
!        | Record indices of the converged wanted Ritz values  |
!        | Mark the select array for possible reordering       |
!        %-----------------------------------------------------%
!
         numcnv = 0
         do 11 j = 1,ncv
            temp1 = max(eps23, abs(workl(irz+ncv-j)) )
            jj = workl(bounds + ncv - j)
            if (numcnv < nconv .and. &
                workl(ibd+jj-1) <= tol*temp1) then
               select(jj) = .true.
               numcnv = numcnv + 1
               if (jj > nev) reord = .true.
            endif
   11    continue
!
!        %-----------------------------------------------------------%
!        | Check the count (numcnv) of converged Ritz values with    |
!        | the number (nconv) reported by _saupd.  If these two      |
!        | are different then there has probably been an error       |
!        | caused by incorrect passing of the _saupd data.           |
!        %-----------------------------------------------------------%
!
         if (msglvl > 2) then
             call ivout(logfil, 1, numcnv, ndigit, &
                  '_seupd: Number of specified eigenvalues')
             call ivout(logfil, 1, nconv, ndigit, &
                  '_seupd: Number of "converged" eigenvalues')
         end if
!
         if (numcnv /= nconv) then
            info = -17
            go to 9000
         end if
!
!        %-----------------------------------------------------------%
!        | Call LAPACK routine _steqr to compute the eigenvalues and |
!        | eigenvectors of the final symmetric tridiagonal matrix H. |
!        | Initialize the eigenvector matrix Q to the identity.      |
!        %-----------------------------------------------------------%
!
         call scopy(ncv-1, workl(ih+1), 1, workl(ihb), 1)
         call scopy(ncv, workl(ih+ldh), 1, workl(ihd), 1)
!
         call ssteqr('Identity', ncv, workl(ihd), workl(ihb), &
                      workl(iq) , ldq, workl(iw), ierr)
!
         if (ierr /= 0) then
            info = -8
            go to 9000
         end if
!
         if (msglvl > 1) then
            call scopy(ncv, workl(iq+ncv-1), ldq, workl(iw), 1)
            call svout(logfil, ncv, workl(ihd), ndigit, &
                '_seupd: NCV Ritz values of the final H matrix')
            call svout(logfil, ncv, workl(iw), ndigit, &
                 '_seupd: last row of the eigenvector matrix for H')
         end if
!
         if (reord) then
!
!           %---------------------------------------------%
!           | Reordered the eigenvalues and eigenvectors  |
!           | computed by _steqr so that the "converged"  |
!           | eigenvalues appear in the first NCONV       |
!           | positions of workl(ihd), and the associated |
!           | eigenvectors appear in the first NCONV      |
!           | columns.                                    |
!           %---------------------------------------------%
!
            leftptr = 1
            rghtptr = ncv
!
            if (ncv == 1) go to 30
!
 20         if (select(leftptr)) then
!
!              %-------------------------------------------%
!              | Search, from the left, for the first Ritz |
!              | value that has not converged.             |
!              %-------------------------------------------%
!
               leftptr = leftptr + 1
!
            else if ( .not. select(rghtptr)) then
!
!              %----------------------------------------------%
!              | Search, from the right, the first Ritz value |
!              | that has converged.                          |
!              %----------------------------------------------%
!
               rghtptr = rghtptr - 1
!
            else
!
!              %----------------------------------------------%
!              | Swap the Ritz value on the left that has not |
!              | converged with the Ritz value on the right   |
!              | that has converged.  Swap the associated     |
!              | eigenvector of the tridiagonal matrix H as   |
!              | well.                                        |
!              %----------------------------------------------%
!
               temp = workl(ihd+leftptr-1)
               workl(ihd+leftptr-1) = workl(ihd+rghtptr-1)
               workl(ihd+rghtptr-1) = temp
               call scopy(ncv, workl(iq+ncv*(leftptr-1)), 1, &
                          workl(iw), 1)
               call scopy(ncv, workl(iq+ncv*(rghtptr-1)), 1, &
                          workl(iq+ncv*(leftptr-1)), 1)
               call scopy(ncv, workl(iw), 1, &
                          workl(iq+ncv*(rghtptr-1)), 1)
               leftptr = leftptr + 1
               rghtptr = rghtptr - 1
!
            end if
!
            if (leftptr < rghtptr) go to 20
!
 30      end if
!
         if (msglvl > 2) then
             call svout (logfil, ncv, workl(ihd), ndigit, &
             '_seupd: The eigenvalues of H--reordered')
         end if
!
!        %----------------------------------------%
!        | Load the converged Ritz values into D. |
!        %----------------------------------------%
!
         call scopy(nconv, workl(ihd), 1, d, 1)
!
      else
!
!        %-----------------------------------------------------%
!        | Ritz vectors not required. Load Ritz values into D. |
!        %-----------------------------------------------------%
!
         call scopy(nconv, workl(ritz), 1, d, 1)
         call scopy(ncv, workl(ritz), 1, workl(ihd), 1)
!
      end if
!
!     %------------------------------------------------------------------%
!     | Transform the Ritz values and possibly vectors and corresponding |
!     | Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values |
!     | (and corresponding data) are returned in ascending order.        |
!     %------------------------------------------------------------------%
!
      if (type == 'REGULR') then
!
!        %---------------------------------------------------------%
!        | Ascending sort of wanted Ritz values, vectors and error |
!        | bounds. Not necessary if only Ritz values are desired.  |
!        %---------------------------------------------------------%
!
         if (rvec) then
            call ssesrt('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call scopy(ncv, workl(bounds), 1, workl(ihb), 1)
         end if
!
      else
!
!        %-------------------------------------------------------------%
!        | *  Make a copy of all the Ritz values.                      |
!        | *  Transform the Ritz values back to the original system.   |
!        |    For TYPE = 'SHIFTI' the transformation is                |
!        |             lambda = 1/theta + sigma                        |
!        |    For TYPE = 'BUCKLE' the transformation is                |
!        |             lambda = sigma * theta / ( theta - 1 )          |
!        |    For TYPE = 'CAYLEY' the transformation is                |
!        |             lambda = sigma * (theta + 1) / (theta - 1 )     |
!        |    where the theta are the Ritz values returned by ssaupd.  |
!        | NOTES:                                                      |
!        | *The Ritz vectors are not affected by the transformation.   |
!        |  They are only reordered.                                   |
!        %-------------------------------------------------------------%
!
         call scopy (ncv, workl(ihd), 1, workl(iw), 1)
         if (type == 'SHIFTI') then
            do 40 k=1, ncv
               workl(ihd+k-1) = one / workl(ihd+k-1) + sigma
  40        continue
         else if (type == 'BUCKLE') then
            do 50 k=1, ncv
               workl(ihd+k-1) = sigma * workl(ihd+k-1) / &
                                (workl(ihd+k-1) - one)
  50        continue
         else if (type == 'CAYLEY') then
            do 60 k=1, ncv
               workl(ihd+k-1) = sigma * (workl(ihd+k-1) + one) / &
                                (workl(ihd+k-1) - one)
  60        continue
         end if
!
!        %-------------------------------------------------------------%
!        | *  Store the wanted NCONV lambda values into D.             |
!        | *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)   |
!        |    into ascending order and apply sort to the NCONV theta   |
!        |    values in the transformed system. We will need this to   |
!        |    compute Ritz estimates in the original system.           |
!        | *  Finally sort the lambda`s into ascending order and apply |
!        |    to Ritz vectors if wanted. Else just sort lambda`s into  |
!        |    ascending order.                                         |
!        | NOTES:                                                      |
!        | *workl(iw:iw+ncv-1) contain the theta ordered so that they  |
!        |  match the ordering of the lambda. We`ll use them again for |
!        |  Ritz vector purification.                                  |
!        %-------------------------------------------------------------%
!
         call scopy(nconv, workl(ihd), 1, d, 1)
         call ssortr('LA', .true., nconv, workl(ihd), workl(iw))
         if (rvec) then
            call ssesrt('LA', rvec , nconv, d, ncv, workl(iq), ldq)
         else
            call scopy(ncv, workl(bounds), 1, workl(ihb), 1)
            call sscal(ncv, bnorm2/rnorm, workl(ihb), 1)
            call ssortr('LA', .true., nconv, d, workl(ihb))
         end if
!
      end if
!
!     %------------------------------------------------%
!     | Compute the Ritz vectors. Transform the wanted |
!     | eigenvectors of the symmetric tridiagonal H by |
!     | the Lanczos basis matrix V.                    |
!     %------------------------------------------------%
!
      if (rvec .and. howmny == 'A') then
!
!        %----------------------------------------------------------%
!        | Compute the QR factorization of the matrix representing  |
!        | the wanted invariant subspace located in the first NCONV |
!        | columns of workl(iq,ldq).                                |
!        %----------------------------------------------------------%
!
         call sgeqr2(ncv, nconv        , workl(iq) , &
                      ldq, workl(iw+ncv), workl(ihb), &
                      ierr)
!
!        %--------------------------------------------------------%
!        | * Postmultiply V by Q.                                 |
!        | * Copy the first NCONV columns of VQ into Z.           |
!        | The N by NCONV matrix Z is now a matrix representation |
!        | of the approximate invariant subspace associated with  |
!        | the Ritz values in workl(ihd).                         |
!        %--------------------------------------------------------%
!
         call sorm2r('Right', 'Notranspose', n        , &
                      ncv    , nconv        , workl(iq), &
                      ldq    , workl(iw+ncv), v        , &
                      ldv    , workd(n+1)   , ierr)
         call slacpy('All', n, nconv, v, ldv, z, ldz)
!
!        %-----------------------------------------------------%
!        | In order to compute the Ritz estimates for the Ritz |
!        | values in both systems, need the last row of the    |
!        | eigenvector matrix. Remember, it`s in factored form |
!        %-----------------------------------------------------%
!
         do 65 j = 1, ncv-1
            workl(ihb+j-1) = zero
  65     continue
         workl(ihb+ncv-1) = one
         call sorm2r('Left', 'Transpose'  , ncv       , &
                      1     , nconv        , workl(iq) , &
                      ldq   , workl(iw+ncv), workl(ihb), &
                      ncv   , temp         , ierr)
!
      else if (rvec .and. howmny == 'S') then
!
!     Not yet implemented. See remark 2 above.
!
      end if
!
      if (type == 'REGULR' .and. rvec) then
!
            do 70 j=1, ncv
               workl(ihb+j-1) = rnorm * abs( workl(ihb+j-1) )
 70         continue
!
      else if (type /= 'REGULR' .and. rvec) then
!
!        %-------------------------------------------------%
!        | *  Determine Ritz estimates of the theta.       |
!        |    If RVEC = .true. then compute Ritz estimates |
!        |               of the theta.                     |
!        |    If RVEC = .false. then copy Ritz estimates   |
!        |              as computed by ssaupd.             |
!        | *  Determine Ritz estimates of the lambda.      |
!        %-------------------------------------------------%
!
         call sscal (ncv, bnorm2, workl(ihb), 1)
         if (type == 'SHIFTI') then

            do 80 k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) ) &
                              / workl(iw+k-1)**2
 80         continue

         else if (type == 'BUCKLE') then

            do 90 k=1, ncv
               workl(ihb+k-1) = sigma * abs( workl(ihb+k-1) ) &
                              / (workl(iw+k-1)-one )**2
 90         continue

         else if (type == 'CAYLEY') then

            do k=1, ncv
               workl(ihb+k-1) = abs( workl(ihb+k-1) &
                              / workl(iw+k-1)*(workl(iw+k-1)-one) )
            end do

         end if

      end if

      if (type /= 'REGULR' .and. msglvl > 1) then
         call svout(logfil, nconv, d, ndigit, &
                '_seupd: Untransformed converged Ritz values')
         call svout(logfil, nconv, workl(ihb), ndigit, &
           '_seupd: Ritz estimates of the untransformed Ritz values')
      else if (msglvl > 1) then
         call svout(logfil, nconv, d, ndigit, &
                '_seupd: Converged Ritz values')
         call svout(logfil, nconv, workl(ihb), ndigit, &
           '_seupd: Associated Ritz estimates')
      end if
!
!     %-------------------------------------------------%
!     | Ritz vector purification step. Formally perform |
!     | one of inverse subspace iteration. Only used    |
!     | for MODE = 3,4,5. See reference 7               |
!     %-------------------------------------------------%
!
      if (rvec .and. (type == 'SHIFTI' .or. type .eq. 'CAYLEY')) then
!
         do 110 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1) &
                        / workl(iw+k)
 110     continue
!
      else if (rvec .and. type == 'BUCKLE') then
!
         do 120 k=0, nconv-1
            workl(iw+k) = workl(iq+k*ldq+ncv-1) &
                        / (workl(iw+k)-one)
 120     continue
!
      end if
!
      if (type /= 'REGULR') &
         call sger (n, nconv, one, resid, 1, workl(iw), 1, z, ldz)
!
 9000 continue
!
      return
!
!     %---------------%
!     | End of sseupd|
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssgets
!
!\Description:
!  Given the eigenvalues of the symmetric tridiagonal matrix H,
!  computes the NP shifts AMU that are zeros of the polynomial of
!  degree NP which filters out components of the unwanted eigenvectors
!  corresponding to the AMU's based on some given criteria.
!
!  NOTE: This is called even in the case of user specified shifts in
!  order to sort the eigenvalues, and error bounds of H for later use.
!
!\Usage:
!  call ssgets
!     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
!
!\Arguments
!  ISHIFT  Integer.  (INPUT)
!          Method for selecting the implicit shifts at each iteration.
!          ISHIFT = 0: user specified shifts
!          ISHIFT = 1: exact shift with respect to the matrix H.
!
!  WHICH   Character*2.  (INPUT)
!          Shift selection criteria.
!          'LM' -> KEV eigenvalues of largest magnitude are retained.
!          'SM' -> KEV eigenvalues of smallest magnitude are retained.
!          'LA' -> KEV eigenvalues of largest value are retained.
!          'SA' -> KEV eigenvalues of smallest value are retained.
!          'BE' -> KEV eigenvalues, half from each end of the spectrum.
!                  If KEV is odd, compute one more from the high end.
!
!  KEV      Integer.  (INPUT)
!          KEV+NP is the size of the matrix H.
!
!  NP      Integer.  (INPUT)
!          Number of implicit shifts to be computed.
!
!  RITZ    Real array of length KEV+NP.  (INPUT/OUTPUT)
!          On INPUT, RITZ contains the eigenvalues of H.
!          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
!          are in the first NP locations and the wanted part is in
!          the last KEV locations.  When exact shifts are selected, the
!          unwanted part corresponds to the shifts to be applied.
!
!  BOUNDS  Real array of length KEV+NP.  (INPUT/OUTPUT)
!          Error bounds corresponding to the ordering in RITZ.
!
!  SHIFTS  Real array of length NP.  (INPUT/OUTPUT)
!          On INPUT:  contains the user specified shifts if ISHIFT = 0.
!          On OUTPUT: contains the shifts sorted into decreasing order
!          of magnitude with respect to the Ritz estimates contained in
!          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     ssortr  ARPACK utility sorting routine.
!     ivout   ARPACK utility routine that prints integers.
!     svout   ARPACK utility routine that prints vectors.
!     scopy   Level 1 BLAS that copies one vector to another.
!     sswap   Level 1 BLAS that swaps the contents of two vectors.
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/93: Version ' 2.1'
!
!\SCCS Information: @(#)
! FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
!
!\Remarks
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssgets ( ishift, which, kev, np, ritz, bounds, shifts )
!
!     %----------------------------------------------------%
!     | Include files for debugging and timing information |
!     %----------------------------------------------------%
!
      include   'debug.h'
      include   'stat.h'
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      integer    ishift, kev, np
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 bounds(kev+np), ritz(kev+np), shifts(np)
!
!     %------------%
!     | Parameters |
!     %------------%
!
      Real &
                 one, zero
      parameter (one = 1.0E+0, zero = 0.0E+0)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    kevd2, msglvl
!
!     %----------------------%
!     | External Subroutines |
!     %----------------------%
!
      external   sswap, scopy, ssortr
!
!     %---------------------%
!     | Intrinsic Functions |
!     %---------------------%
!
      intrinsic    max, min
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
!     %-------------------------------%
!     | Initialize timing statistics  |
!     | & message level for debugging |
!     %-------------------------------%
!
      call cpu_time (t0)
      msglvl = msgets
!
      if (which == 'BE') then
!
!        %-----------------------------------------------------%
!        | Both ends of the spectrum are requested.            |
!        | Sort the eigenvalues into algebraically increasing  |
!        | order first then swap high end of the spectrum next |
!        | to low end in appropriate locations.                |
!        | NOTE: when np < floor(kev/2) be careful not to swap |
!        | overlapping locations.                              |
!        %-----------------------------------------------------%
!
         call ssortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2
         if ( kev > 1 ) then
            call sswap ( min(kevd2,np), ritz, 1, &
                         ritz( max(kevd2,np)+1 ), 1)
            call sswap ( min(kevd2,np), bounds, 1, &
                         bounds( max(kevd2,np)+1 ), 1)
         end if
!
      else
!
!        %----------------------------------------------------%
!        | LM, SM, LA, SA case.                               |
!        | Sort the eigenvalues of H into the desired order   |
!        | and apply the resulting order to BOUNDS.           |
!        | The eigenvalues are sorted so that the wanted part |
!        | are always in the last KEV locations.               |
!        %----------------------------------------------------%
!
         call ssortr (which, .true., kev+np, ritz, bounds)
      end if
!
      if (ishift == 1 .and. np > 0) then
!
!        %-------------------------------------------------------%
!        | Sort the unwanted Ritz values used as shifts so that  |
!        | the ones with largest Ritz estimates are first.       |
!        | This will tend to minimize the effects of the         |
!        | forward instability of the iteration when the shifts  |
!        | are applied in subroutine ssapps.                     |
!        %-------------------------------------------------------%
!
         call ssortr ('SM', .true., np, bounds, ritz)
         call scopy (np, ritz, 1, shifts, 1)
      end if
!
      call cpu_time (t1)
      tsgets = tsgets + (t1 - t0)
!
      if (msglvl > 0) then
         call ivout (logfil, 1, kev, ndigit, '_sgets: KEV is')
         call ivout (logfil, 1, np, ndigit, '_sgets: NP is')
         call svout (logfil, kev+np, ritz, ndigit, &
              '_sgets: Eigenvalues of current H matrix')
         call svout (logfil, kev+np, bounds, ndigit, &
              '_sgets: Associated Ritz estimates')
      end if
!
      return
!
!     %---------------%
!     | End of ssgets |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssortc
!
!\Description:
!  Sorts the complex array in XREAL and XIMAG into the order
!  specified by WHICH and optionally applies the permutation to the
!  real array Y. It is assumed that if an element of XIMAG is
!  nonzero, then its negative is also an element. In other words,
!  both members of a complex conjugate pair are to be sorted and the
!  pairs are kept adjacent to each other.
!
!\Usage:
!  call ssortc
!     ( WHICH, APPLY, N, XREAL, XIMAG, Y )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> sort XREAL,XIMAG into increasing order of magnitude.
!          'SM' -> sort XREAL,XIMAG into decreasing order of magnitude.
!          'LR' -> sort XREAL into increasing order of algebraic.
!          'SR' -> sort XREAL into decreasing order of algebraic.
!          'LI' -> sort XIMAG into increasing order of magnitude.
!          'SI' -> sort XIMAG into decreasing order of magnitude.
!          NOTE: If an element of XIMAG is non-zero, then its negative
!                is also an element.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to array Y.
!          APPLY = .FALSE. -> do not apply the sorted order to array Y.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  XREAL,  Real array of length N.  (INPUT/OUTPUT)
!  XIMAG   Real and imaginary part of the array to be sorted.
!
!  Y       Real array of length N.  (INPUT/OUTPUT)
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     xx/xx/92: Version ' 2.1'
!               Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#)
! FILE: sortc.F   SID: 2.3   DATE OF SID: 4/20/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssortc (which, apply, n, xreal, ximag, y)
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      logical    apply
      integer    n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 xreal(0:n-1), ximag(0:n-1), y(0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, igap, j
      Real &
                 temp, temp1, temp2
!
!     %--------------------%
!     | External Functions |
!     %--------------------%
!
      Real &
                 slapy2
      external   slapy2
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      igap = n / 2
!
      if (which == 'LM') then
!
!        %------------------------------------------------------%
!        | Sort XREAL,XIMAG into increasing order of magnitude. |
!        %------------------------------------------------------%
!
   10    continue
         if (igap == 0) go to 9000
!
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if (j<0) go to 30
!
            temp1 = slapy2(xreal(j),ximag(j))
            temp2 = slapy2(xreal(j+igap),ximag(j+igap))
!
            if (temp1>temp2) then
                temp = xreal(j)
                xreal(j) = xreal(j+igap)
                xreal(j+igap) = temp
!
                temp = ximag(j)
                ximag(j) = ximag(j+igap)
                ximag(j+igap) = temp
!
                if (apply) then
                    temp = y(j)
                    y(j) = y(j+igap)
                    y(j+igap) = temp
                end if
            else
                go to 30
            end if
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
!
      else if (which == 'SM') then
!
!        %------------------------------------------------------%
!        | Sort XREAL,XIMAG into decreasing order of magnitude. |
!        %------------------------------------------------------%
!
   40    continue
         if (igap == 0) go to 9000
!
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j < 0) go to 60
!
            temp1 = slapy2(xreal(j),ximag(j))
            temp2 = slapy2(xreal(j+igap),ximag(j+igap))
!
            if (temp1<temp2) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
!
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
!
      else if (which == 'LR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into increasing order of algebraic. |
!        %------------------------------------------------%
!
   70    continue
         if (igap == 0) go to 9000
!
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j<0) go to 90
!
            if (xreal(j)>xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
!
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
!
      else if (which == 'SR') then
!
!        %------------------------------------------------%
!        | Sort XREAL into decreasing order of algebraic. |
!        %------------------------------------------------%
!
  100    continue
         if (igap == 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j<0) go to 120
!
            if (xreal(j)<xreal(j+igap)) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
!
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
!
      else if (which == 'LI') then
!
!        %------------------------------------------------%
!        | Sort XIMAG into increasing order of magnitude. |
!        %------------------------------------------------%
!
  130    continue
         if (igap == 0) go to 9000
         do 150 i = igap, n-1
            j = i-igap
  140       continue
!
            if (j<0) go to 150
!
            if (abs(ximag(j))>abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
!
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 150
            endif
            j = j-igap
            go to 140
  150    continue
         igap = igap / 2
         go to 130
!
      else if (which == 'SI') then
!
!        %------------------------------------------------%
!        | Sort XIMAG into decreasing order of magnitude. |
!        %------------------------------------------------%
!
  160    continue
         if (igap == 0) go to 9000
         do 180 i = igap, n-1
            j = i-igap
  170       continue
!
            if (j<0) go to 180
!
            if (abs(ximag(j))<abs(ximag(j+igap))) then
               temp = xreal(j)
               xreal(j) = xreal(j+igap)
               xreal(j+igap) = temp
!
               temp = ximag(j)
               ximag(j) = ximag(j+igap)
               ximag(j+igap) = temp
!
               if (apply) then
                  temp = y(j)
                  y(j) = y(j+igap)
                  y(j+igap) = temp
               end if
            else
               go to 180
            endif
            j = j-igap
            go to 170
  180    continue
         igap = igap / 2
         go to 160
      end if
!
 9000 continue
      return
!
!     %---------------%
!     | End of ssortc |
!     %---------------%
!
      end
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: ssortr
!
!\Description:
!  Sort the array X1 in the order specified by WHICH and optionally
!  applies the permutation to the array X2.
!
!\Usage:
!  call ssortr
!     ( WHICH, APPLY, N, X1, X2 )
!
!\Arguments
!  WHICH   Character*2.  (Input)
!          'LM' -> X1 is sorted into increasing order of magnitude.
!          'SM' -> X1 is sorted into decreasing order of magnitude.
!          'LA' -> X1 is sorted into increasing order of algebraic.
!          'SA' -> X1 is sorted into decreasing order of algebraic.
!
!  APPLY   Logical.  (Input)
!          APPLY = .TRUE.  -> apply the sorted order to X2.
!          APPLY = .FALSE. -> do not apply the sorted order to X2.
!
!  N       Integer.  (INPUT)
!          Size of the arrays.
!
!  X1      Real array of length N.  (INPUT/OUTPUT)
!          The array to be sorted.
!
!  X2      Real array of length N.  (INPUT/OUTPUT)
!          Only referenced if APPLY = .TRUE.
!
!\EndDoc
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Author
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\Revision history:
!     12/16/93: Version ' 2.1'.
!               Adapted from the sort routine in LANSO.
!
!\SCCS Information: @(#)
! FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
!
!\EndLib
!
!-----------------------------------------------------------------------
!
      subroutine ssortr (which, apply, n, x1, x2)
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      character*2 which
      logical    apply
      integer    n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 x1(0:n-1), x2(0:n-1)
!
!     %---------------%
!     | Local Scalars |
!     %---------------%
!
      integer    i, igap, j
      Real &
                 temp
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      igap = n / 2
!
      if (which == 'SA') then
!
!        X1 is sorted into decreasing order of algebraic.
!
   10    continue
         if (igap == 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
!
            if (j<0) go to 30
!
            if (x1(j)<x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
!
      else if (which == 'SM') then
!
!        X1 is sorted into decreasing order of magnitude.
!
   40    continue
         if (igap == 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
!
            if (j<0) go to 60
!
            if (abs(x1(j))<abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
!
      else if (which == 'LA') then
!
!        X1 is sorted into increasing order of algebraic.
!
   70    continue
         if (igap == 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
!
            if (j<0) go to 90
!
            if (x1(j)>x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
!
      else if (which == 'LM') then
!
!        X1 is sorted into increasing order of magnitude.
!
  100    continue
         if (igap == 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
!
            if (j<0) go to 120
!
            if (abs(x1(j))>abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
!
 9000 continue
      return
!
!     %---------------%
!     | End of ssortr |
!     %---------------%
!
      end
      subroutine sstatn
!
!
!! SSTATN initializes the statistics for the nonsymmetric Arnolid iteration.
!
!     %--------------------------------%
!     | See stat.doc for documentation |
!     %--------------------------------%
!
      include   'stat.h'
!
!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%
!
      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0
!
      tnaupd = 0.0E+0
      tnaup2 = 0.0E+0
      tnaitr = 0.0E+0
      tneigh = 0.0E+0
      tngets = 0.0E+0
      tnapps = 0.0E+0
      tnconv = 0.0E+0
      titref = 0.0E+0
      tgetv0 = 0.0E+0
      trvec  = 0.0E+0
!
!     %----------------------------------------------------%
!     | User time including reverse communication overhead |
!     %----------------------------------------------------%
!
      tmvopx = 0.0E+0
      tmvbx  = 0.0E+0
!
      return
!
!
!     %---------------%
!     | End of sstatn |
!     %---------------%
!
      end
      subroutine sstats
!
!! SSTATS initializes the statistics for the symmetric Arnolid iteration.
!
!     %--------------------------------%
!     | See stat.doc for documentation |
!     %--------------------------------%
      include   'stat.h'

!     %-----------------------%
!     | Executable Statements |
!     %-----------------------%

      nopx   = 0
      nbx    = 0
      nrorth = 0
      nitref = 0
      nrstrt = 0

      tsaupd = 0.0E+0
      tsaup2 = 0.0E+0
      tsaitr = 0.0E+0
      tseigt = 0.0E+0
      tsgets = 0.0E+0
      tsapps = 0.0E+0
      tsconv = 0.0E+0
      titref = 0.0E+0
      tgetv0 = 0.0E+0
      trvec  = 0.0E+0

!     %----------------------------------------------------%
!     | User time including reverse communication overhead |
!     %----------------------------------------------------%
      tmvopx = 0.0E+0
      tmvbx  = 0.0E+0

      return
!
!     End of sstats
!
      end
      subroutine sstqrb ( n, d, e, z, work, info )
!
!! SSTQRB computes all eigenvalues of a symmetric tridiagonal matrix.
!
!-----------------------------------------------------------------------
!\BeginDoc
!
!\Name: sstqrb
!
!\Description:
!  Computes all eigenvalues and the last component of the eigenvectors
!  of a symmetric tridiagonal matrix using the implicit QL or QR method.
!
!  This is mostly a modification of the LAPACK routine ssteqr.
!  See Remarks.
!
!\Usage:
!  call sstqrb
!     ( N, D, E, Z, WORK, INFO )
!
!\Arguments
!  N       Integer.  (INPUT)
!          The number of rows and columns in the matrix.  N >= 0.
!
!  D       Real array, dimension (N).  (INPUT/OUTPUT)
!          On entry, D contains the diagonal elements of the
!          tridiagonal matrix.
!          On exit, D contains the eigenvalues, in ascending order.
!          If an error exit is made, the eigenvalues are correct
!          for indices 1,2,...,INFO-1, but they are unordered and
!          may not be the smallest eigenvalues of the matrix.
!
!  E       Real array, dimension (N-1).  (INPUT/OUTPUT)
!          On entry, E contains the subdiagonal elements of the
!          tridiagonal matrix in positions 1 through N-1.
!          On exit, E has been destroyed.
!
!  Z       Real array, dimension (N).  (OUTPUT)
!          On exit, Z contains the last row of the orthonormal
!          eigenvector matrix of the symmetric tridiagonal matrix.
!          If an error exit is made, Z contains the last row of the
!          eigenvector matrix associated with the stored eigenvalues.
!
!  WORK    Real array, dimension (max(1,2*N-2)).  (WORKSPACE)
!          Workspace used in accumulating the transformation for
!          computing the last components of the eigenvectors.
!
!  INFO    Integer.  (OUTPUT)
!          = 0:  normal return.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = +i, the i-th eigenvalue has not converged
!                              after a total of  30*N  iterations.
!
!\Remarks
!  1. None.
!
!-----------------------------------------------------------------------
!
!\BeginLib
!
!\Local variables:
!     xxxxxx  real
!
!\Routines called:
!     saxpy   Level 1 BLAS that computes a vector triad.
!     scopy   Level 1 BLAS that copies one vector to another.
!     sswap   Level 1 BLAS that swaps the contents of two vectors.
!     lsame   LAPACK character comparison routine.
!     slae2   LAPACK routine that computes the eigenvalues of a 2-by-2
!             symmetric matrix.
!     slaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric
!             matrix.
!     slamch  LAPACK routine that determines machine constants.
!     slanst  LAPACK routine that computes the norm of a matrix.
!     slapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
!     slartg  LAPACK Givens rotation construction routine.
!     slascl  LAPACK routine for careful scaling of a matrix.
!     slaset  LAPACK matrix initialization routine.
!     slasr   LAPACK routine that applies an orthogonal transformation to
!             a matrix.
!     slasrt  LAPACK sorting routine.
!     ssteqr  LAPACK routine that computes eigenvalues and eigenvectors
!             of a symmetric tridiagonal matrix.
!     xerbla  LAPACK error handler routine.
!
!\Authors
!     Danny Sorensen               Phuong Vu
!     Richard Lehoucq              CRPC / Rice University
!     Dept. of Computational &     Houston, Texas
!     Applied Mathematics
!     Rice University
!     Houston, Texas
!
!\SCCS Information: @(#)
! FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
!
!\Remarks
!     1. Starting with version 2.5, this routine is a modified version
!        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
!        only commeted out and new lines inserted.
!        All lines commented out have "c$$$" at the beginning.
!        Note that the LAPACK version 1.0 subroutine SSTEQR contained
!        bugs.
!
!\EndLib
!
!-----------------------------------------------------------------------
!
!
!     %------------------%
!     | Scalar Arguments |
!     %------------------%
!
      integer    info, n
!
!     %-----------------%
!     | Array Arguments |
!     %-----------------%
!
      Real &
                 d( n ), e( n-1 ), z( n ), work( 2*n-2 )
!
!     .. parameters ..
      Real &
                         zero, one, two, three
      parameter          ( zero = 0.0E+0, one = 1.0E+0, &
                           two = 2.0E+0, three = 3.0E+0 )
      integer            maxit
      parameter          ( maxit = 30 )
!     ..
!     .. local scalars ..
      integer            i, icompz, ii, iscale, j, jtot, k, l, l1, lend, &
                         lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1, &
                         nm1, nmaxit
      Real &
                         anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2, &
                         s, safmax, safmin, ssfmax, ssfmin, tst
!     ..
!     .. external functions ..
      logical            lsame
      Real &
                         slamch, slanst, slapy2
      external           lsame, slamch, slanst, slapy2
!     ..
!     .. external subroutines ..
      external           slae2, slaev2, slartg, slascl, slaset, slasr, &
                         slasrt, sswap, xerbla
!     ..
!     .. intrinsic functions ..
      intrinsic          abs, max, sign, sqrt
!     ..
!     .. executable statements ..
!
!     test the input parameters.
!
      info = 0
!
!$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
!$$$         ICOMPZ = 0
!$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
!$$$         ICOMPZ = 1
!$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
!$$$         ICOMPZ = 2
!$$$      ELSE
!$$$         ICOMPZ = -1
!$$$      END IF
!$$$      IF( ICOMPZ.LT.0 ) THEN
!$$$         INFO = -1
!$$$      ELSE IF( N.LT.0 ) THEN
!$$$         INFO = -2
!$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
!$$$     $         N ) ) ) THEN
!$$$         INFO = -6
!$$$      END IF
!$$$      IF( INFO.NE.0 ) THEN
!$$$         CALL XERBLA( 'SSTEQR', -INFO )
!$$$         RETURN
!$$$      END IF
!
!    *** New starting with version 2.5 ***
!
      icompz = 2
!    *************************************
!
!     quick return if possible
!
      if( n==0 ) &
         return
!
      if( n==1 ) then
         if( icompz==2 )  z( 1 ) = one
         return
      end if
!
!     determine the unit roundoff and over/underflow thresholds.
!
      eps = slamch( 'e' )
      eps2 = eps**2
      safmin = slamch( 's' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2
!
!     compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
!$$      if( icompz==2 )
!$$$     $   call slaset( 'full', n, n, zero, one, z, ldz )
!
!     *** New starting with version 2.5 ***
!
      if ( icompz == 2 ) then
         do 5 j = 1, n-1
            z(j) = zero
  5      continue
         z( n ) = one
      end if
!     *************************************
!
      nmaxit = n*maxit
      jtot = 0
!
!     determine where the matrix splits and choose ql or qr iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      l1 = 1
      nm1 = n - 1
!
   10 continue
      if( l1>n ) &
         go to 160
      if( l1>1 ) &
         e( l1-1 ) = zero
      if( l1<=nm1 ) then
         do m = l1, nm1
            tst = abs( e( m ) )
            if( tst==zero ) &
               go to 30
            if( tst<=( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+ &
                1 ) ) ) )*eps ) then
               e( m ) = zero
               go to 30
            end if
         end do
      end if
      m = n
!
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if( lend==l ) &
         go to 10
!
!     scale submatrix in rows and columns l to lend
!
      anorm = slanst( 'i', lend-l+1, d( l ), e( l ) )
      iscale = 0
      if( anorm==zero ) &
         go to 10
      if( anorm>ssfmax ) then
         iscale = 1
         call slascl( 'g', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n, &
                      info )
         call slascl( 'g', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n, &
                      info )
      else if( anorm<ssfmin ) then
         iscale = 2
         call slascl( 'g', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n, &
                      info )
         call slascl( 'g', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n, &
                      info )
      end if
!
!     choose between ql and qr iteration
!
      if( abs( d( lend ) )<abs( d( l ) ) ) then
         lend = lsv
         l = lendsv
      end if
!
      if( lend>l ) then
!
!        ql iteration
!
!        look for small subdiagonal element.
!
   40    continue
         if( l/=lend ) then
            lendm1 = lend - 1
            do 50 m = l, lendm1
               tst = abs( e( m ) )**2
               if( tst<=( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+ &
                   safmin )go to 60
   50       continue
         end if
!
         m = lend
!
   60    continue
         if( m<lend ) &
            e( m ) = zero
         p = d( l )
         if( m==l ) &
            go to 80
!
!        if remaining matrix is 2-by-2, use slae2 or slaev2
!        to compute its eigensystem.
!
         if( m==l+1 ) then
            if( icompz>0 ) then
               call slaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
               work( l ) = c
               work( n-1+l ) = s
!$$$               call slasr( 'r', 'v', 'b', n, 2, work( l ),
!$$$     $                     work( n-1+l ), z( 1, l ), ldz )
!
!              *** New starting with version 2.5 ***
!
               tst      = z(l+1)
               z(l+1) = c*tst - s*z(l)
               z(l)   = s*tst + c*z(l)
!              *************************************
            else
               call slae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
            end if
            d( l ) = rt1
            d( l+1 ) = rt2
            e( l ) = zero
            l = l + 2
            if( l<=lend ) &
               go to 40
            go to 140
         end if
!
         if( jtot==nmaxit ) &
            go to 140
         jtot = jtot + 1
!
!        form shift.
!
         g = ( d( l+1 )-p ) / ( two*e( l ) )
         r = slapy2( g, one )
         g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
!
         s = one
         c = one
         p = zero
!
!        inner loop
!
         mm1 = m - 1
         do 70 i = mm1, l, -1
            f = s*e( i )
            b = c*e( i )
            call slartg( g, f, c, s, r )
            if( i/=m-1 ) &
               e( i+1 ) = r
            g = d( i+1 ) - p
            r = ( d( i )-g )*s + two*c*b
            p = s*r
            d( i+1 ) = g + p
            g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
            if( icompz>0 ) then
               work( i ) = c
               work( n-1+i ) = -s
            end if
!
   70    continue
!
!        if eigenvectors are desired, then apply saved rotations.
!
         if( icompz>0 ) then
            mm = m - l + 1
!$$$            call slasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ),
!$$$     $                  z( 1, l ), ldz )
!
!             *** New starting with version 2.5 ***
!
              call slasr( 'r', 'v', 'b', 1, mm, work( l ), &
                          work( n-1+l ), z( l ), 1 )
!             *************************************
         end if
!
         d( l ) = d( l ) - p
         e( l ) = g
         go to 40
!
!        eigenvalue found.
!
   80    continue
         d( l ) = p
!
         l = l + 1
         if( l<=lend ) &
            go to 40
         go to 140
!
      else
!
!        qr iteration
!
!        look for small superdiagonal element.
!
   90    continue
         if( l/=lend ) then
            lendp1 = lend + 1
            do m = l, lendp1, -1
               tst = abs( e( m-1 ) )**2
               if( tst<=( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+ &
                   safmin )go to 110
            end do
         end if

         m = lend

  110    continue
         if( m>lend ) &
            e( m-1 ) = zero
         p = d( l )
         if( m==l ) &
            go to 130
!
!        if remaining matrix is 2-by-2, use slae2 or slaev2
!        to compute its eigensystem.
!
         if( m==l-1 ) then
            if( icompz>0 ) then
               call slaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
!$$$               work( m ) = c
!$$$               work( n-1+m ) = s
!$$$               call slasr( 'r', 'v', 'f', n, 2, work( m ),
!$$$     $                     work( n-1+m ), z( 1, l-1 ), ldz )
!
!               *** New starting with version 2.5 ***
!
                tst      = z(l)
                z(l)   = c*tst - s*z(l-1)
                z(l-1) = s*tst + c*z(l-1)
!               *************************************
            else
               call slae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
            end if
            d( l-1 ) = rt1
            d( l ) = rt2
            e( l-1 ) = zero
            l = l - 2
            if( l>=lend ) &
               go to 90
            go to 140
         end if
!
         if( jtot==nmaxit ) &
            go to 140
         jtot = jtot + 1
!
!        form shift.
!
         g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
         r = slapy2( g, one )
         g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )
!
         s = one
         c = one
         p = zero
!
!        inner loop
!
         lm1 = l - 1
         do 120 i = m, lm1
            f = s*e( i )
            b = c*e( i )
            call slartg( g, f, c, s, r )
            if( i/=m ) &
               e( i-1 ) = r
            g = d( i ) - p
            r = ( d( i+1 )-g )*s + two*c*b
            p = s*r
            d( i ) = g + p
            g = c*r - b
!
!           if eigenvectors are desired, then save rotations.
!
            if( icompz>0 ) then
               work( i ) = c
               work( n-1+i ) = s
            end if
!
  120    continue
!
!        if eigenvectors are desired, then apply saved rotations.
!
         if( icompz>0 ) then
            mm = l - m + 1
!$$$            call slasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ),
!$$$     $                  z( 1, m ), ldz )
!
!           *** New starting with version 2.5 ***
!
            call slasr( 'r', 'v', 'f', 1, mm, work( m ), work( n-1+m ), &
                        z( m ), 1 )
!           *************************************
         end if
!
         d( l ) = d( l ) - p
         e( lm1 ) = g
         go to 90
!
!        eigenvalue found.
!
  130    continue
         d( l ) = p
!
         l = l - 1
         if( l>=lend ) &
            go to 90
         go to 140
!
      end if
!
!     undo scaling if necessary
!
  140 continue
      if( iscale==1 ) then
         call slascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1, &
                      d( lsv ), n, info )
         call slascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ), &
                      n, info )
      else if( iscale==2 ) then
         call slascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1, &
                      d( lsv ), n, info )
         call slascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ), &
                      n, info )
      end if
!
!     check for no convergence to an eigenvalue after a total
!     of n*maxit iterations.
!
      if( jtot<nmaxit ) &
         go to 10
      do 150 i = 1, n - 1
         if( e( i )/=zero ) &
            info = info + 1
  150 continue
      go to 190
!
!     order eigenvalues and eigenvectors.
!
  160 continue
      if( icompz==0 ) then
!
!        use quick sort
!
         call slasrt( 'i', n, d, info )
!
      else
!
!        use selection sort to minimize swaps of eigenvectors
!
         do 180 ii = 2, n
            i = ii - 1
            k = i
            p = d( i )
            do 170 j = ii, n
               if( d( j )<p ) then
                  k = j
                  p = d( j )
               end if
  170       continue
            if( k/=i ) then
               d( k ) = d( i )
               d( i ) = p
!$$$               call sswap( n, z( 1, i ), 1, z( 1, k ), 1 )
!           *** New starting with version 2.5 ***
!
               p    = z(k)
               z(k) = z(i)
               z(i) = p
!           *************************************
            end if
  180    continue
      end if
!
  190 continue
      return
!
!     %---------------%
!     | End of sstqrb |
!     %---------------%
!
      end
      integer function icnteq (n, array, value)
!
!! ICNTEQ counts the number of vector elements equal to a given value.
!
      integer    n, value
      integer    array(*)
!
      k = 0
      do i = 1, n
         if (array(i) == value) k = k + 1
      end do
      icnteq = k
!
      return
      end
      subroutine icopy( n, lx, incx, ly, incy )
!--------------------------------------------------------------------
!
!! ICOPY copies an integer vector lx to an integer vector ly.
!
!\Usage:
!     call icopy ( n, lx, inc, ly, incy )
!
!\Arguments:
!    n        integer (input)
!             On entry, n is the number of elements of lx to be
!             copied to ly.
!
!    lx       integer array (input)
!             On entry, lx is the integer vector to be copied.
!
!    incx     integer (input)
!             On entry, incx is the increment between elements of lx.
!
!    ly       integer array (input)
!             On exit, ly is the integer vector that contains the
!             copy of lx.
!
!    incy     integer (input)
!             On entry, incy is the increment between elements of ly.
!
!\Enddoc
!
!--------------------------------------------------------------------
!
!     ----------------------------
!     Specifications for arguments
!     ----------------------------
      integer    incx, incy, n
      integer    lx( 1 ), ly( 1 )
!
!     ----------------------------------
!     Specifications for local variables
!     ----------------------------------
      integer           i, ix, iy
!
!     --------------------------
!     First executable statement
!     --------------------------
      if( n<=0 ) &
         return
      if( incx==1 .and. incy.eq.1 ) &
         go to 20
!
!.....code for unequal increments or equal increments
!     not equal to 1
      ix = 1
      iy = 1
      if( incx<0 ) &
         ix = ( -n+1 )*incx + 1
      if( incy<0 ) &
         iy = ( -n+1 )*incy + 1
      do i = 1, n
         ly( iy ) = lx( ix )
         ix = ix + incx
         iy = iy + incy
      end do

      return
!
!.....code for both increments equal to 1
!
   20 continue
      do 30 i = 1, n
         ly( i ) = lx( i )
   30 continue
      return
      end
      subroutine iset (n, value, array, inc)
!
!! ISET sets all entries of an integer vector to a scalar.
!
      integer    n, value, inc
      integer    array(*)
!
      array(1:n) = value
!
      return
      end
      subroutine iswap (n,sx,incx,sy,incy)
!
!! ISWAP interchanges two vectors.
!     uses unrolled loops for increments equal to 1.
!     jack dongarra, linpack, 3/11/78.
!
      integer sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n<=0)return
      if(incx==1.and.incy.eq.1)go to 20
!
!       code for unequal increments or equal increments not equal
!         to 1
!
      ix = 1
      iy = 1
      if(incx<0)ix = (-n+1)*incx + 1
      if(incy<0)iy = (-n+1)*incy + 1
      do i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
      end do

      return
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
   20 m = mod(n,3)
      if( m == 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n < 3 ) return
   40 mp1 = m + 1

      do i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
      end do

      return
      end
      SUBROUTINE IVOUT (LOUT, N, IX, IDIGIT, IFMT)
!-----------------------------------------------------------------------
!
!! IVOUT is an integer vector output routine.
!
!  Usage:      CALL IVOUT (LOUT, N, IX, IDIGIT, IFMT)
!
!  Arguments
!     N      - Length of array IX. (Input)
!     IX     - Integer array to be printed. (Input)
!     IFMT   - Format to be used in printing array IX. (Input)
!     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input)
!              If IDIGIT .LT. 0, printing is done with 72 columns.
!              If IDIGIT .GT. 0, printing is done with 132 columns.
!
!-----------------------------------------------------------------------
!
!     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IX(*), N, IDIGIT, LOUT
      CHARACTER  IFMT*(*)
!     ...
!     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, NDIGIT, K1, K2, LLL
      CHARACTER*80 LINE
!     ...
!     ... SPECIFICATIONS INTRINSICS
      INTRINSIC          MIN
!
!
      LLL = MIN ( LEN ( IFMT ), 80 )
      DO I = 1, LLL
          LINE(I:I) = '-'
      end do

      DO I = LLL+1, 80
          LINE(I:I) = ' '
      end do

      WRITE ( LOUT, 2000 ) IFMT, LINE(1:LLL)
 2000 FORMAT ( /1X, A  /1X, A )

      IF (N .LE. 0) RETURN
      NDIGIT = IDIGIT
      IF (IDIGIT .EQ. 0) NDIGIT = 4
!
!=======================================================================
!             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
!=======================================================================
!
      IF (IDIGIT .LT. 0) THEN
!
      NDIGIT = -IDIGIT
      IF (NDIGIT .LE. 4) THEN
         DO K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
         end do

      ELSE IF (NDIGIT .LE. 6) THEN
         DO K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
         end do

      ELSE IF (NDIGIT .LE. 10) THEN
         DO 50 K1 = 1, N, 5
            K2 = MIN0(N,K1+4)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
   50    CONTINUE

      ELSE
         DO 70 K1 = 1, N, 3
            K2 = MIN0(N,K1+2)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
   70    CONTINUE
      END IF
!
!=======================================================================
!             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
!=======================================================================
!
      ELSE

      IF (NDIGIT .LE. 4) THEN
         DO 90 K1 = 1, N, 20
            K2 = MIN0(N,K1+19)
            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
   90    CONTINUE

      ELSE IF (NDIGIT .LE. 6) THEN
         DO K1 = 1, N, 15
            K2 = MIN0(N,K1+14)
            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
         end do

      ELSE IF (NDIGIT .LE. 10) THEN
         DO 130 K1 = 1, N, 10
            K2 = MIN0(N,K1+9)
            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
  130    CONTINUE

      ELSE
         DO 150 K1 = 1, N, 7
            K2 = MIN0(N,K1+6)
            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
  150    CONTINUE
      END IF
      END IF
      WRITE (LOUT,1004)

 1000 FORMAT(1X,I4,' - ',I4,':',20(1X,I5))
 1001 FORMAT(1X,I4,' - ',I4,':',15(1X,I7))
 1002 FORMAT(1X,I4,' - ',I4,':',10(1X,I11))
 1003 FORMAT(1X,I4,' - ',I4,':',7(1X,I15))
 1004 FORMAT(1X,' ')

      RETURN
      END
      SUBROUTINE SMOUT( LOUT, M, N, A, LDA, IDIGIT, IFMT )
!
!-----------------------------------------------------------------------
!
!! SMOUT is a real matrix output routine.
!
!  Usage:      CALL SMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT)
!
!  Arguments
!     M      - Number of rows of A.  (Input)
!     N      - Number of columns of A.  (Input)
!     A      - Real M by N matrix to be printed.  (Input)
!     LDA    - Leading dimension of A exactly as specified in the
!              dimension statement of the calling program.  (Input)
!     IFMT   - Format to be used in printing matrix A.  (Input)
!     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
!              If IDIGIT .LT. 0, printing is done with 72 columns.
!              If IDIGIT .GT. 0, printing is done with 132 columns.
!
!-----------------------------------------------------------------------
!
!     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER            M, N, IDIGIT, LDA, LOUT
      REAL               A( LDA, * )
      CHARACTER          IFMT*( * )
!     ...
!     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I, J, NDIGIT, K1, K2, LLL
      CHARACTER*1        ICOL( 3 )
      CHARACTER*80       LINE
!     ...
!     ... SPECIFICATIONS INTRINSICS
      INTRINSIC          MIN
!
      DATA               ICOL( 1 ), ICOL( 2 ), ICOL( 3 ) / 'C', 'o', &
                         'l' /
!
      LLL = MIN( LEN( IFMT ), 80 )
      DO I = 1, LLL
         LINE( I: I ) = '-'
      end do

      DO i = LLL + 1, 80
         LINE( I: I ) = ' '
      end do

      WRITE( LOUT, '(/1x,a/1x,a)' )IFMT, LINE( 1: LLL )

      IF( M.LE.0 .OR. N.LE.0 .OR. LDA.LE.0 ) &
         RETURN
      NDIGIT = IDIGIT
      IF( IDIGIT.EQ.0 ) &
         NDIGIT = 4
!
!=======================================================================
!             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
!=======================================================================
!
      IF( IDIGIT.LT.0 ) THEN
         NDIGIT = -IDIGIT
         IF( NDIGIT.LE.4 ) THEN
            DO 40 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, 9998 )( ICOL, I, I = K1, K2 )
               DO 30 I = 1, M
                  WRITE( LOUT, 9994 )I, ( A( I, J ), J = K1, K2 )
   30          CONTINUE
   40       CONTINUE
!
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 60 K1 = 1, N, 4
               K2 = MIN0( N, K1+3 )
               WRITE( LOUT, 9997 )( ICOL, I, I = K1, K2 )
               DO 50 I = 1, M
                  WRITE( LOUT, 9993 )I, ( A( I, J ), J = K1, K2 )
   50          CONTINUE
   60       CONTINUE
!
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 80 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
               WRITE( LOUT, 9996 )( ICOL, I, I = K1, K2 )
               DO 70 I = 1, M
                  WRITE( LOUT, 9992 )I, ( A( I, J ), J = K1, K2 )
   70          CONTINUE
   80       CONTINUE

         ELSE
            DO K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, 9995 )( ICOL, I, I = K1, K2 )
               DO 90 I = 1, M
                  WRITE( LOUT, 9991 )I, ( A( I, J ), J = K1, K2 )
   90          CONTINUE
            end do
         END IF
!
!=======================================================================
!             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
!=======================================================================
!
      ELSE
         IF( NDIGIT.LE.4 ) THEN
            DO 120 K1 = 1, N, 10
               K2 = MIN0( N, K1+9 )
               WRITE( LOUT, 9998 )( ICOL, I, I = K1, K2 )
               DO 110 I = 1, M
                  WRITE( LOUT, 9994 )I, ( A( I, J ), J = K1, K2 )
  110          CONTINUE
  120       CONTINUE
!
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 140 K1 = 1, N, 8
               K2 = MIN0( N, K1+7 )
               WRITE( LOUT, 9997 )( ICOL, I, I = K1, K2 )
               DO 130 I = 1, M
                  WRITE( LOUT, 9993 )I, ( A( I, J ), J = K1, K2 )
  130          CONTINUE
  140       CONTINUE
!
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 160 K1 = 1, N, 6
               K2 = MIN0( N, K1+5 )
               WRITE( LOUT, 9996 )( ICOL, I, I = K1, K2 )
               DO 150 I = 1, M
                  WRITE( LOUT, 9992 )I, ( A( I, J ), J = K1, K2 )
  150          CONTINUE
  160       CONTINUE
!
         ELSE
            DO 180 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, 9995 )( ICOL, I, I = K1, K2 )
               DO 170 I = 1, M
                  WRITE( LOUT, 9991 )I, ( A( I, J ), J = K1, K2 )
  170          CONTINUE
  180       CONTINUE
         END IF
      END IF
      WRITE( LOUT, 9990 )
!
 9998 FORMAT( 10X, 10( 4X, 3A1, I4, 1X ) )
 9997 FORMAT( 10X, 8( 5X, 3A1, I4, 2X ) )
 9996 FORMAT( 10X, 6( 7X, 3A1, I4, 4X ) )
 9995 FORMAT( 10X, 5( 9X, 3A1, I4, 6X ) )
 9994 FORMAT( 1X, ' Row', I4, ':', 1X, 1P10E12.3 )
 9993 FORMAT( 1X, ' Row', I4, ':', 1X, 1P8E14.5 )
 9992 FORMAT( 1X, ' Row', I4, ':', 1X, 1P6E18.9 )
 9991 FORMAT( 1X, ' Row', I4, ':', 1X, 1P5E22.13 )
 9990 FORMAT( 1X, ' ' )
!
      RETURN
      END
      SUBROUTINE SVOUT( LOUT, N, SX, IDIGIT, IFMT )
!-----------------------------------------------------------------------
!
!! SVOUT is a real vector output routine.
!
!  Usage:      CALL SVOUT (LOUT, N, SX, IDIGIT, IFMT)
!
!  Arguments
!     N      - Length of array SX.  (Input)
!     SX     - Real array to be printed.  (Input)
!     IFMT   - Format to be used in printing array SX.  (Input)
!     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
!              If IDIGIT .LT. 0, printing is done with 72 columns.
!              If IDIGIT .GT. 0, printing is done with 132 columns.
!
!     ...
!     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER            N, IDIGIT, LOUT
      REAL               SX( * )
      CHARACTER          IFMT*( * )
!     ...
!     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER            I, NDIGIT, K1, K2, LLL
      CHARACTER*80       LINE
!     ...
!     ... FIRST EXECUTABLE STATEMENT
!
!
      LLL = MIN( LEN( IFMT ), 80 )
      DO I = 1, LLL
         LINE( I: I ) = '-'
      end do

      DO 20 I = LLL + 1, 80
         LINE( I: I ) = ' '
   20 CONTINUE
!
      WRITE( LOUT, '(/1x,a/1x,a)' )IFMT, LINE( 1: LLL )


      IF( N.LE.0 ) &
         RETURN
      NDIGIT = IDIGIT
      IF( IDIGIT.EQ.0 ) &
         NDIGIT = 4
!
!=======================================================================
!             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
!=======================================================================
!
      IF( IDIGIT.LT.0 ) THEN
         NDIGIT = -IDIGIT
         IF( NDIGIT.LE.4 ) THEN
            DO 30 K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, 9998 )K1, K2, ( SX( I ), I = K1, K2 )
   30       CONTINUE
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 40 K1 = 1, N, 4
               K2 = MIN0( N, K1+3 )
               WRITE( LOUT, 9997 )K1, K2, ( SX( I ), I = K1, K2 )
   40       CONTINUE
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 50 K1 = 1, N, 3
               K2 = MIN0( N, K1+2 )
               WRITE( LOUT, 9996 )K1, K2, ( SX( I ), I = K1, K2 )
   50       CONTINUE
         ELSE
            DO 60 K1 = 1, N, 2
               K2 = MIN0( N, K1+1 )
               WRITE( LOUT, 9995 )K1, K2, ( SX( I ), I = K1, K2 )
   60       CONTINUE
         END IF
!
!=======================================================================
!             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
!=======================================================================
!
      ELSE
         IF( NDIGIT.LE.4 ) THEN
            DO 70 K1 = 1, N, 10
               K2 = MIN0( N, K1+9 )
               WRITE( LOUT, 9998 )K1, K2, ( SX( I ), I = K1, K2 )
   70       CONTINUE
         ELSE IF( NDIGIT.LE.6 ) THEN
            DO 80 K1 = 1, N, 8
               K2 = MIN0( N, K1+7 )
               WRITE( LOUT, 9997 )K1, K2, ( SX( I ), I = K1, K2 )
   80       CONTINUE
         ELSE IF( NDIGIT.LE.10 ) THEN
            DO 90 K1 = 1, N, 6
               K2 = MIN0( N, K1+5 )
               WRITE( LOUT, 9996 )K1, K2, ( SX( I ), I = K1, K2 )
   90       CONTINUE
         ELSE
            DO K1 = 1, N, 5
               K2 = MIN0( N, K1+4 )
               WRITE( LOUT, 9995 )K1, K2, ( SX( I ), I = K1, K2 )
            end do
         END IF
      END IF
      WRITE( LOUT, 9994 )
      RETURN
 9998 FORMAT( 1X, I4, ' - ', I4, ':', 1P10E12.3 )
 9997 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P8E14.5 )
 9996 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P6E18.9 )
 9995 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P5E24.13 )
 9994 FORMAT( 1X, ' ' )
      END
subroutine timestamp ( )
!
!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
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
!
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
