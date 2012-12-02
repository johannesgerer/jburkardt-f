program main

!*****************************************************************************80
!
!! MAIN is the main program for PBLAS_PRB.
!
!  Discussion:
!
!    PBLAS_PRB sets matrix descriptors and calls the PBLAS routines.
!
!
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!
!  Author:
!
!    Antoine Petitet, August 1995 (petitet@cs.utk.edu)
!
  implicit none
!
  integer, parameter :: dblesz = 8
  integer, parameter :: totmem = 200000


 integer            memsiz
  parameter          ( memsiz = totmem / dblesz )
  integer            block_cyclic_2d, csrc_, ctxt_, dlen_, dt_, &
                     lld_, mb_, m_, nb_, n_, rsrc_
  parameter          ( block_cyclic_2d = 1, dlen_ = 9, dt_ = 1, &
                       ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6, &
                       rsrc_ = 7, csrc_ = 8, lld_ = 9 )

  integer            iam, iaseed, ibseed, icseed, ictxt, info, ipa, &
                     ipb, ipc, ipw, k, kp, kq, m, mp, my_col, my_row, &
                     n, nb, nout, npcol, nprocs, nprow, nq, worksiz

  double precision bnrm2
  integer desca(dlen_)
  integer descb(dlen_)
  integer descc(dlen_ )
  double precision mem(memsiz)
  double precision, parameter :: one = 1.0d+00
  character ( len = 80 ) outfile

  external           blacs_exit, blacs_get, blacs_gridexit, &
                     blacs_gridinfo, blacs_gridinit, blacs_pinfo, &
                     descinit, igsum2d, pdmatgen, pdpblasinfo, &
                     pdnrm2, pdgemv, pdgemm, pdlaprnt

  integer            numroc
  external           numroc 
!
!  Get starting information.
!
  call BLACS_PINFO ( iam, NPROCS )

  call PDPBLASINFO ( OUTFILE, nout, M, N, K, nb, NPROW, NPCOL, mem, &
                    iam, NPROCS )
!
!  Define process grid.
!
  call BLACS_GET( -1, 0, ICTXT )
  call BLACS_GRIDINIT( ICTXT, 'Row-major', NPROW, NPCOL )
  call BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, my_row, my_col )
!
!  Go to bottom of process grid loop if this case doesn't use my
!  process.
!
  if ( NPROW <= my_row .or. NPCOL <= my_col ) then
     go to 20
  end if

  MP = numroc( M, nb, my_row, 0, NPROW )
  KP = numroc( K, nb, my_row, 0, NPROW )
  KQ = numroc( K, nb, my_col, 0, NPCOL )
  nq = numroc( N, nb, my_col, 0, NPCOL )
!
!  Initialize the array descriptor for the matrix A, B and C.
!
  call descinit( desca, M, K, nb, nb, 0, 0, ICTXT, MAX( 1, MP ), &
                 INFO )
  call descinit( descb, K, N, nb, nb, 0, 0, ICTXT, MAX( 1, KP ), &
                 INFO )
  call descinit( descc, M, N, nb, nb, 0, 0, ICTXT, MAX( 1, MP ), &
                 INFO )
!
!  Assign pointers into mem for SCALAPACK arrays, A is
!  allocated starting at position mem( 1 ).
!
  IPA = 1
  IPB = IPA + desca( LLD_ )*KQ
  IPC = IPB + descb( LLD_ )*NQ
  IPW = IPC + descc( LLD_ )*NQ

  WORKSIZ = nb
!
!  Check for adequate memory for problem size.
!
  INFO = 0
  if ( IPW+WORKSIZ > memSIZ ) then
     if ( iam == 0 ) then
        write ( nout, FMT = 9998 ) 'test', ( IPW+WORKSIZ )*DBLESZ
     end if
     INFO = 1
  end if
!
!  Check all processes for an error.
!
  call IGSUM2D( ICTXT, 'All', ' ', 1, 1, INFO, 1, -1, 0 )

  if ( INFO > 0 ) then
    if ( iam == 0 ) then
      write ( nout, FMT = 9999 ) 'memORY'
    end if
    go to 10
  end if
!
!  Generate random matrices A, B and C.
!
  IASEED = 100
  call PDMATGEN( ICTXT, 'No transpose', 'No transpose', desca( M_ ), &
                 desca( N_ ), DESCA( MB_ ), DESCA( nb_ ), &
                 mem( IPA ), desca( LLD_ ), DESCA( RSRC_ ), &
                 desca( CSRC_ ), IASEED, 0, MP, 0, KQ, my_row, my_col, &
                 NPROW, NPCOL )

  IBSEED = 200
  call PDMATGEN( ICTXT, 'No transpose', 'No transpose', descb( M_ ), &
                 descb( N_ ), DESCB( MB_ ), DESCB( nb_ ), &
                 mem( IPB ), descb( LLD_ ), DESCB( RSRC_ ), &
                 descb( CSRC_ ), IBSEED, 0, KP, 0, NQ, my_row, my_col, &
                 NPROW, NPCOL )

  ICSEED = 300
  call PDMATGEN( ICTXT, 'No transpose', 'No transpose', descc( M_ ), &
                 descc( N_ ), DESCC( MB_ ), DESCC( nb_ ), &
                 mem( IPC ), descc( LLD_ ), DESCC( RSRC_ ), &
                 descc( CSRC_ ), ICSEED, 0, MP, 0, NQ, my_row, my_col, &
                 NPROW, NPCOL )
!
!*********************************************************************
!     Call Level 1 PBLAS routine
!*********************************************************************
!
  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) &
           '***********************************************'
     write ( nout, FMT = * ) &
           'Example of Level 1 PBLAS routine call: (PDNRM2)'
     write ( nout, FMT = * ) &
           '***********************************************'
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Vector X:'
     write ( nout, FMT = * )
  end if

  call PDLAPRNT( K, 1, mem( IPB ), 1, 1, descb, 0, 0, &
                 'X', nout, mem( IPW ) )

  call PDNRM2( K, BNRM2, mem( IPB ), 1, 1, descb, 1 )

  if ( my_row == 0 .and. my_col == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) '|| X || = sqrt( X''*X ) = ', BNRM2
     write ( nout, FMT = * )
  end if
!
!*********************************************************************
!     Call Level 2 PBLAS routine
!*********************************************************************
!
  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) &
           '***********************************************'
     write ( nout, FMT = * ) &
           'Example of Level 2 PBLAS routine call: (PDGEMV)'
     write ( nout, FMT = * ) &
           '***********************************************'
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Matrix A:'
     write ( nout, FMT = * )
  end if

  call PDLAPRNT( M, K, mem( IPA ), 1, 1, desca, 0, 0, &
                 'A', nout, mem( IPW ) )

  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Vector X:'
     write ( nout, FMT = * )
  end if
  call PDLAPRNT( K, 1, mem( IPB ), 1, 1, descb, 0, 0, &
                 'X', nout, mem( IPW ) )

  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Vector Y:'
     write ( nout, FMT = * )
  end if
  call PDLAPRNT( M, 1, mem( IPC ), 1, 1, descc, 0, 0, &
                 'Y', nout, mem( IPW ) )

  call PDGEMV( 'No transpose', M, K, ONE, mem( IPA ), 1, 1, desca, &
               mem( IPB ), 1, 1, descb, 1, ONE, MEM( IPC ), 1, 1, &
               descc, 1 )

  if ( my_row == 0 .and. my_col == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Y := Y + A * X'
     write ( nout, FMT = * )
  end if

  call PDLAPRNT( M, 1, mem( IPC ), 1, 1, descc, 0, 0, &
                 'Y', nout, mem( IPW ) )
!
!*********************************************************************
!     Call Level 3 PBLAS routine
!*********************************************************************
!
  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) &
           '***********************************************'
     write ( nout, FMT = * ) &
           'Example of Level 3 PBLAS routine call: (PDGEMM)'
     write ( nout, FMT = * ) &
           '***********************************************'
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Matrix A:'
     write ( nout, FMT = * )
  end if

  call PDLAPRNT( M, K, mem( IPA ), 1, 1, desca, 0, 0, &
                 'A', nout, mem( IPW ) )

  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Matrix B:'
     write ( nout, FMT = * )
  end if

  call PDLAPRNT( K, N, mem( IPB ), 1, 1, descb, 0, 0, &
                 'B', nout, mem( IPW ) )

  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' Matrix C:'
     write ( nout, FMT = * )
  end if
  call PDLAPRNT( M, N, mem( IPC ), 1, 1, descc, 0, 0, &
                 'C', nout, mem( IPW ) )

  call PDGEMM( 'No transpose', 'No transpose', M, N, K, ONE, &
               mem( IPA ), 1, 1, desca, MEM( IPB ), 1, 1, descb, &
               ONE, mem( IPC ), 1, 1, descc )

  if ( my_row == 0 .and. my_col == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * ) ' C := C + A * B'
     write ( nout, FMT = * )
  end if

  call PDLAPRNT( M, N, mem( IPC ), 1, 1, descc, 0, 0, &
                 'C', nout, mem( IPW ) )

   10 CONTINUE

  call BLACS_GRIDEXIT ( ICTXT )

   20 CONTINUE
!
!  Print ending messages and close output file.
!
  if ( iam == 0 ) then
     write ( nout, FMT = * )
     write ( nout, FMT = * )
     write ( nout, FMT = 9997 )
     write ( nout, FMT = * )
     if ( nout /= 6 .and. nout /= 0 ) then
        CLOSE ( nout )
     end if
  end if

  call BLACS_EXIT ( 0 )

 9999 FORMAT( 'Bad ', A6, ' parameters: going on to next test case.' )
 9998 FORMAT( 'Unable to perform ', A, ': need TOTMEM of at least', &
          I11 )
 9997 FORMAT( 'END OF TESTS.' )

  stop
end
SUBROUTINE PDMATGEN( ICTXT, AFORM, DIAG, M, N, MB, NB, A, LDA, &
                       IAROW, IACOL, ISEED, IROFF, IRNUM, ICOFF, &
                       ICNUM, MYROW, MYCOL, NPROW, NPCOL )
!
!*****************************************************************************80
!
!  Purpose
!  =======
!
!  PDMATGEN : Parallel Real Double precision MATrix GENerator.
!  Generate (or regenerate) a distributed matrix A (or sub-matrix of A).
!
!  Arguments
!  =========
!
!  ICTXT   (global input) INTEGER
!          The BLACS context handle, indicating the global context of
!          the operation. The context itself is global.
!
!  AFORM   (global input) CHARACTER
!          if AFORM = 'S' : A is returned is a symmetric matrix.
!          if AFORM = 'H' : A is returned is a Hermitian matrix.
!          if AFORM = 'T' : A is overwritten with the transpose of
!                           what would normally be generated.
!          if AFORM = 'C' : A is overwritten with the conjugate trans-
!                           pose of what would normally be generated.
!          otherwise a random matrix is generated.
!
!  DIAG    (global input) CHARACTER
!          if DIAG = 'D' : A is diagonally dominant.
!
!  M       (global input) INTEGER
!          The number of rows in the generated distributed matrix.
!
!  N       (global input) INTEGER
!          The number of columns in the generated distributed
!          matrix.
!
!  MB      (global input) INTEGER
!          The row blocking factor of the distributed matrix A.
!
!  NB      (global input) INTEGER
!          The column blocking factor of the distributed matrix A.
!
!  A       (local output) DOUBLE PRECISION, pointer into the local
!          memory to an array of dimension ( LDA, * ) containing the
!          local pieces of the distributed matrix.
!
!  LDA     (local input) INTEGER
!          The leading dimension of the array containing the local
!          pieces of the distributed matrix A.
!
!  IAROW   (global input) INTEGER
!          The row processor coordinate which holds the first block
!          of the distributed matrix A.
!
!  IACOL   (global input) INTEGER
!          The column processor coordinate which holds the first
!          block of the distributed matrix A.
!
!  ISEED   (global input) INTEGER
!          The seed number to generate the distributed matrix A.
!
!  IROFF   (local input) INTEGER
!          The number of local rows of A that have already been
!          generated.  It should be a multiple of MB.
!
!  IRNUM   (local input) INTEGER
!          The number of local rows to be generated.
!
!  ICOFF   (local input) INTEGER
!          The number of local columns of A that have already been
!          generated.  It should be a multiple of NB.
!
!  ICNUM   (local input) INTEGER
!          The number of local columns to be generated.
!
!  MYROW   (local input) INTEGER
!          The row process coordinate of the calling process.
!
!  MYCOL   (local input) INTEGER
!          The column process coordinate of the calling process.
!
!  NPROW   (global input) INTEGER
!          The number of process rows in the grid.
!
!  NPCOL   (global input) INTEGER
!          The number of process columns in the grid.
!
!  Notes
!  =====
!
!  The code is originally developed by David Walker, ORNL,
!  and modified by Jaeyoung Choi, ORNL.
!
!  Reference: G. Fox et al.
!  Section 12.3 of "Solving problems on concurrent processors Vol. I"
!
  CHARACTER AFORM
  character DIAG
  INTEGER            IACOL, IAROW, ICNUM, ICOFF, ICTXT, IRNUM, &
                     IROFF, ISEED, LDA, M, MB, MYCOL, MYROW, N, &
                     NB, NPCOL, NPROW
  DOUBLE PRECISION   A( LDA, * )
  INTEGER            MULT0, MULT1, IADD0, IADD1
  PARAMETER        ( MULT0=20077, MULT1=16838, IADD0=12345, &
                     IADD1=0 )
  DOUBLE PRECISION   ONE, TWO
  PARAMETER          ( ONE = 1.0D+0, TWO = 2.0D+0 )
!     ..
!     .. Local Scalars ..
  LOGICAL            SYMM, HERM, TRAN
  INTEGER            I, IC, IK, INFO, IOFFC, IOFFR, IR, J, JK, &
                     JUMP1, JUMP2, JUMP3, JUMP4, JUMP5, JUMP6, &
                     JUMP7, MAXMN, MEND, MOFF, MP, MRCOL, MRROW, &
                     NEND, NOFF, NPMB, NQ, NQNB
!     ..
!     .. Local Arrays ..
  INTEGER            IADD(2), IA1(2), IA2(2), IA3(2), IA4(2), &
                     IA5(2), IB1(2), IB2(2), IB3(2), IC1(2), IC2(2), &
                     IC3(2), IC4(2), IC5(2), IRAN1(2), IRAN2(2), &
                     IRAN3(2), IRAN4(2), ITMP1(2), ITMP2(2), &
                     ITMP3(2), JSEED(2), MULT(2)
!     ..
!     .. External Subroutines ..
  EXTERNAL           JUMPIT, PXERBLA, SETRAN, XJUMPM
!     ..
!     .. External Functions ..
  LOGICAL            LSAME
  INTEGER            ICEIL, NUMROC
  DOUBLE PRECISION   PDRAND
  EXTERNAL           ICEIL, NUMROC, LSAME, PDRAND
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
  MP   = NUMROC( M, MB, MYROW, IAROW, NPROW )
  NQ   = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
  SYMM = LSAME( AFORM, 'S' )
  HERM = LSAME( AFORM, 'H' )
  TRAN = LSAME( AFORM, 'T' )
!
  INFO = 0
  if ( .NOT.( SYMM.OR.HERM.OR.TRAN ) .AND. &
      .NOT.LSAME( AFORM, 'C' ) .AND. &
      .NOT.LSAME( AFORM, 'N' )            ) THEN
     INFO = 2
  else if ( .NOT.LSAME( DIAG, 'D' ) .AND. &
           .NOT.LSAME( DIAG, 'N' )        ) THEN
     INFO = 3
  else if ( SYMM.OR.HERM ) THEN
     if ( M.NE.N ) THEN
        INFO = 5
     else if ( MB.NE.NB ) THEN
        INFO = 7
     end if
  else if ( M.LT.0 ) THEN
     INFO = 4
  else if ( N.LT.0 ) THEN
     INFO = 5
  else if ( MB.LT.1 ) THEN
     INFO = 6
  else if ( NB.LT.1 ) THEN
     INFO = 7
  else if ( LDA.LT.0 ) THEN
     INFO = 9
  else if ( ( IAROW.LT.0 ).OR.( IAROW.GE.NPROW ) ) THEN
     INFO = 10
  else if ( ( IACOL.LT.0 ).OR.( IACOL.GE.NPCOL ) ) THEN
     INFO = 11
  else if ( MOD(IROFF,MB) > 0 ) THEN
     INFO = 13
  else if ( IRNUM > (MP-IROFF) ) THEN
     INFO = 14
  else if ( MOD(ICOFF,NB) > 0 ) THEN
     INFO = 15
  else if ( ICNUM > (NQ-ICOFF) ) THEN
     INFO = 16
  else if ( ( MYROW.LT.0 ).OR.( MYROW.GE.NPROW ) ) THEN
     INFO = 17
  else if ( ( MYCOL.LT.0 ).OR.( MYCOL.GE.NPCOL ) ) THEN
     INFO = 18
  end if
  if ( INFO.NE.0 ) THEN
     CALL PXERBLA( ICTXT, 'PDMATGEN', INFO )
     RETURN
  end if
!
  MRROW = MOD( NPROW+MYROW-IAROW, NPROW )
  MRCOL = MOD( NPCOL+MYCOL-IACOL, NPCOL )
  NPMB  = NPROW * MB
  NQNB  = NPCOL * NB
  MOFF  = IROFF / MB
  NOFF  = ICOFF / NB
  MEND  = ICEIL(IRNUM, MB) + MOFF
  NEND  = ICEIL(ICNUM, NB) + NOFF
!
  MULT(1)  = MULT0
  MULT(2)  = MULT1
  IADD(1)  = IADD0
  IADD(2)  = IADD1
  JSEED(1) = ISEED
  JSEED(2) = 0
!
!  Symmetric or Hermitian matrix will be generated.
!
  if ( SYMM .OR. HERM ) THEN
!
!        First, generate the lower triangular part (with diagonal block)
!
     JUMP1 = 1
     JUMP2 = NPMB
     JUMP3 = M
     JUMP4 = NQNB
     JUMP5 = NB
     JUMP6 = MRCOL
     JUMP7 = MB*MRROW
!
     CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
     CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
     CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
     CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
     CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
     CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
     CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
     CALL XJUMPM( NOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
     CALL XJUMPM( MOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
     CALL SETRAN( IRAN1, IA1,  IC1 )

     DO I = 1, 2
        IB1(I) = IRAN1(I)
        IB2(I) = IRAN1(I)
        IB3(I) = IRAN1(I)
     end do

     JK = 1
     DO 80 IC = NOFF+1, NEND
        IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
        DO 70 I = 1, NB
           if ( JK  >  ICNUM ) GO TO 90

           IK = 1
           DO 50 IR = MOFF+1, MEND
              IOFFR = ((IR-1)*NPROW+MRROW) * MB

              if ( IOFFR  >  IOFFC ) THEN
                 DO J = 1, MB
                    if ( IK  >  IRNUM ) GO TO 60
                    A(IK,JK) = ONE - TWO*PDRAND(0)
                    IK = IK + 1
                 end do

              else if ( IOFFC .EQ. IOFFR ) THEN
                 IK = IK + I - 1
                 if ( IK  >  IRNUM ) GO TO 60
                 DO J = 1, I-1
                    A(IK,JK) = ONE - TWO*PDRAND(0)
                 end do
                 A(IK,JK) = ONE - TWO*PDRAND(0)
                 DO 40 J = 1, MB-I
                    if ( IK+J  >  IRNUM ) GO TO 60
                    A(IK+J,JK) = ONE - TWO*PDRAND(0)
                    A(IK,JK+J) = A(IK+J,JK)
   40                CONTINUE
                 IK = IK + MB - I + 1
              else
                 IK = IK + MB
              end if
!
              CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
              IB1(1) = IRAN2(1)
              IB1(2) = IRAN2(2)
   50          CONTINUE
!
   60          CONTINUE
           JK = JK + 1
           CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
           IB1(1) = IRAN3(1)
           IB1(2) = IRAN3(2)
           IB2(1) = IRAN3(1)
           IB2(2) = IRAN3(2)
   70       CONTINUE
!
        CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
        IB1(1) = IRAN4(1)
        IB1(2) = IRAN4(2)
        IB2(1) = IRAN4(1)
        IB2(2) = IRAN4(2)
        IB3(1) = IRAN4(1)
        IB3(2) = IRAN4(2)
   80    CONTINUE
!
!  Next, generate the upper triangular part.
!
   90    CONTINUE
     MULT(1)  = MULT0
     MULT(2)  = MULT1
     IADD(1)  = IADD0
     IADD(2)  = IADD1
     JSEED(1) = ISEED
     JSEED(2) = 0

     JUMP1 = 1
     JUMP2 = NQNB
     JUMP3 = N
     JUMP4 = NPMB
     JUMP5 = MB
     JUMP6 = MRROW
     JUMP7 = NB*MRCOL

     CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
     CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
     CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
     CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
     CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
     CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
     CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
     CALL XJUMPM( MOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
     CALL XJUMPM( NOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
     CALL SETRAN( IRAN1, IA1,  IC1 )
!
     DO 100 I = 1, 2
        IB1(I) = IRAN1(I)
        IB2(I) = IRAN1(I)
        IB3(I) = IRAN1(I)
  100    CONTINUE
!
     IK = 1
     DO 150 IR = MOFF+1, MEND
        IOFFR = ((IR-1)*NPROW+MRROW) * MB
        DO 140 J = 1, MB
           if ( IK  >  IRNUM ) GO TO 160
           JK = 1
           DO 120 IC = NOFF+1, NEND
              IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
              if ( IOFFC  >  IOFFR ) THEN
                 DO 110 I = 1, NB
                    if ( JK  >  ICNUM ) GO TO 130
                    A(IK,JK) = ONE - TWO*PDRAND(0)
                    JK = JK + 1
  110                CONTINUE
              else
                 JK = JK + NB
              end if
              CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
              IB1(1) = IRAN2(1)
              IB1(2) = IRAN2(2)
  120          CONTINUE
!
  130          CONTINUE
           IK = IK + 1
           CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
           IB1(1) = IRAN3(1)
           IB1(2) = IRAN3(2)
           IB2(1) = IRAN3(1)
           IB2(2) = IRAN3(2)
  140       CONTINUE
!
        CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
        IB1(1) = IRAN4(1)
        IB1(2) = IRAN4(2)
        IB2(1) = IRAN4(1)
        IB2(2) = IRAN4(2)
        IB3(1) = IRAN4(1)
        IB3(2) = IRAN4(2)
  150    CONTINUE
  160    CONTINUE
!
!     (Conjugate) Transposed matrix A will be generated.
!
  else if ( TRAN .or. LSAME( AFORM, 'C' ) ) THEN
!
     JUMP1 = 1
     JUMP2 = NQNB
     JUMP3 = N
     JUMP4 = NPMB
     JUMP5 = MB
     JUMP6 = MRROW
     JUMP7 = NB*MRCOL
!
     CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
     CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
     CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
     CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
     CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
     CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
     CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
     CALL XJUMPM( MOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
     CALL XJUMPM( NOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
     CALL SETRAN( IRAN1, IA1,  IC1 )
!
     DO 170 I = 1, 2
        IB1(I) = IRAN1(I)
        IB2(I) = IRAN1(I)
        IB3(I) = IRAN1(I)
  170    CONTINUE
!
     IK = 1
     DO 220 IR = MOFF+1, MEND
        IOFFR = ((IR-1)*NPROW+MRROW) * MB
        DO 210 J = 1, MB
           if ( IK  >  IRNUM ) GO TO 230
           JK = 1
           DO 190 IC = NOFF+1, NEND
              IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
              DO 180 I = 1, NB
                 if ( JK  >  ICNUM ) GO TO 200
                 A(IK,JK) = ONE - TWO*PDRAND(0)
                 JK = JK + 1
  180             CONTINUE
              CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
              IB1(1) = IRAN2(1)
              IB1(2) = IRAN2(2)
  190          CONTINUE
!
  200          CONTINUE
           IK = IK + 1
           CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
           IB1(1) = IRAN3(1)
           IB1(2) = IRAN3(2)
           IB2(1) = IRAN3(1)
           IB2(2) = IRAN3(2)
  210       CONTINUE
!
        CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
        IB1(1) = IRAN4(1)
        IB1(2) = IRAN4(2)
        IB2(1) = IRAN4(1)
        IB2(2) = IRAN4(2)
        IB3(1) = IRAN4(1)
        IB3(2) = IRAN4(2)
  220    CONTINUE
  230    CONTINUE
!
!     A random matrix is generated.
!
  else
!
     JUMP1 = 1
     JUMP2 = NPMB
     JUMP3 = M
     JUMP4 = NQNB
     JUMP5 = NB
     JUMP6 = MRCOL
     JUMP7 = MB*MRROW

     CALL XJUMPM( JUMP1, MULT, IADD, JSEED, IRAN1, IA1,   IC1 )
     CALL XJUMPM( JUMP2, MULT, IADD, IRAN1, ITMP1, IA2,   IC2 )
     CALL XJUMPM( JUMP3, MULT, IADD, IRAN1, ITMP1, IA3,   IC3 )
     CALL XJUMPM( JUMP4, IA3,  IC3,  IRAN1, ITMP1, IA4,   IC4 )
     CALL XJUMPM( JUMP5, IA3,  IC3,  IRAN1, ITMP1, IA5,   IC5 )
     CALL XJUMPM( JUMP6, IA5,  IC5,  IRAN1, ITMP3, ITMP1, ITMP2 )
     CALL XJUMPM( JUMP7, MULT, IADD, ITMP3, IRAN1, ITMP1, ITMP2 )
     CALL XJUMPM( NOFF,  IA4,  IC4,  IRAN1, ITMP1, ITMP2, ITMP3 )
     CALL XJUMPM( MOFF,  IA2,  IC2,  ITMP1, IRAN1, ITMP2, ITMP3 )
     CALL SETRAN( IRAN1, IA1,  IC1 )

     DO I = 1, 2
        IB1(I) = IRAN1(I)
        IB2(I) = IRAN1(I)
        IB3(I) = IRAN1(I)
     end do

     JK = 1
     DO 290 IC = NOFF+1, NEND
        IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
        DO 280 I = 1, NB
           if ( JK  >  ICNUM ) GO TO 300
           IK = 1
           DO 260 IR = MOFF+1, MEND
              IOFFR = ((IR-1)*NPROW+MRROW) * MB
              DO 250 J = 1, MB
                 if ( IK  >  IRNUM ) GO TO 270
                 A(IK,JK) = ONE - TWO*PDRAND(0)
                 IK = IK + 1
  250             CONTINUE
              CALL JUMPIT( IA2, IC2, IB1, IRAN2 )
              IB1(1) = IRAN2(1)
              IB1(2) = IRAN2(2)
  260          CONTINUE

  270          CONTINUE
           JK = JK + 1
           CALL JUMPIT( IA3, IC3, IB2, IRAN3 )
           IB1(1) = IRAN3(1)
           IB1(2) = IRAN3(2)
           IB2(1) = IRAN3(1)
           IB2(2) = IRAN3(2)
  280       CONTINUE

        CALL JUMPIT( IA4, IC4, IB3, IRAN4 )
        IB1(1) = IRAN4(1)
        IB1(2) = IRAN4(2)
        IB2(1) = IRAN4(1)
        IB2(2) = IRAN4(2)
        IB3(1) = IRAN4(1)
        IB3(2) = IRAN4(2)
  290    CONTINUE
  300    CONTINUE
  end if
!
!     Diagonally dominant matrix will be generated.
!
  if ( LSAME( DIAG, 'D' ) ) THEN
     if ( MB.NE.NB ) THEN
        WRITE(*,*) 'Diagonally dominant matrices with rowNB not'// &
                   ' equal colNB is not supported!'
        RETURN
     end if
!
     MAXMN = MAX(M, N)
     JK    = 1
     DO 340 IC = NOFF+1, NEND
        IOFFC = ((IC-1)*NPCOL+MRCOL) * NB
        IK    = 1
        DO 320 IR = MOFF+1, MEND
           IOFFR = ((IR-1)*NPROW+MRROW) * MB
           if ( IOFFC.EQ.IOFFR ) THEN
              DO J = 0, MB-1
                 if ( IK  >  IRNUM ) GO TO 330
                 A(IK,JK+J) = ABS(A(IK,JK+J)) + MAXMN
                 IK = IK + 1
              end do
           else
              IK = IK + MB
           end if
  320       CONTINUE
  330       CONTINUE
        JK = JK + NB
  340    CONTINUE
  end if

  return
end
SUBROUTINE PDPBLASINFO ( SUMMRY, NOUT, M, N, K, NB, NPROW, NPCOL, &
  WORK, IAM, NPROCS )
!
!*****************************************************************************80
!
!  -- PBLAS example code --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!
!     Written by Antoine Petitet, August 1995 (petitet@cs.utk.edu)
!
!     This program shows how to set the matrix descriptors and call
!     the PBLAS routines.
!
!     .. Scalar Arguments ..
  CHARACTER ( len = * ) SUMMRY
  INTEGER            IAM, K, M, N, NB, NOUT, NPCOL, NPROCS, NPROW
!     ..
!     .. Array Arguments ..
  INTEGER WORK( * )
!     ..
!
!     .. Parameters ..
  integer, PARAMETER :: NIN = 11
!     ..
  CHARACTER ( len = 79 )      USRINFO
  INTEGER            ICTXT
!     ..
!     .. External Subroutines ..
  EXTERNAL           BLACS_ABORT, BLACS_GET, BLACS_GRIDEXIT, &
                     BLACS_GRIDINIT, BLACS_SETUP, IGEBR2D, IGEBS2D
!     ..
!     .. Executable Statements ..
!
!     Process 0 reads the input data, broadcasts to other processes and
!     writes needed information to NOUT
!
  if ( IAM == 0 ) THEN
!
!        Open file and skip data file header
!
     OPEN( NIN, FILE='pblas.dat', STATUS='OLD' )
     READ( NIN, FMT = * ) SUMMRY
     SUMMRY = ' '
!
!        Read in user-supplied info about machine type, compiler, etc.
!
     READ( NIN, FMT = 9999 ) USRINFO
!
!        Read name and unit number for summary output file
!
     READ( NIN, FMT = * ) SUMMRY
     READ( NIN, FMT = * ) NOUT
     if ( NOUT /= 0 .AND. NOUT /= 6 ) &
        OPEN( NOUT, FILE = SUMMRY, STATUS = 'UNKNOWN' )
!
!        Read and check the parameter values for the tests.
!
!        Get matrix dimensions
!
     READ( NIN, FMT = * ) M
     READ( NIN, FMT = * ) N
     READ( NIN, FMT = * ) K
!
!  Get value of NB
!
     READ( NIN, FMT = * ) NB
!
!        Get grid shape
!
     READ( NIN, FMT = * ) NPROW
     READ( NIN, FMT = * ) NPCOL
!
!        Close input file
!
     CLOSE( NIN )
!
!        If underlying system needs additional set up, do it now
!
     if ( NPROCS < 1 ) THEN
        NPROCS = NPROW * NPCOL
        CALL BLACS_SETUP( IAM, NPROCS )
     end if
!
!        Temporarily define blacs grid to include all processes so
!        information can be broadcast to all processes
!
     CALL BLACS_GET( -1, 0, ICTXT )
     CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )
!
!        Pack information arrays and broadcast
!
     WORK( 1 ) = M
     WORK( 2 ) = N
     WORK( 3 ) = K
     WORK( 4 ) = NB
     WORK( 5 ) = NPROW
     WORK( 6 ) = NPCOL
     CALL IGEBS2D( ICTXT, 'All', ' ', 6, 1, WORK, 6 )
!
!        regurgitate input
!
     WRITE( NOUT, FMT = 9999 ) &
                 'PBLAS Examples driver.'
     WRITE( NOUT, FMT = 9999 ) USRINFO
     WRITE( NOUT, FMT = * )
     WRITE( NOUT, FMT = 9999 ) &
                 'The matrices A, B and C are randomly '// &
                 'generated for each test.'
     WRITE( NOUT, FMT = * )
     WRITE( NOUT, FMT = 9999 ) &
                 'An explanation of the input/output '// &
                 'parameters follows:'
!
     WRITE( NOUT, FMT = 9999 ) &
                 'M       : The number of rows in the '// &
                 'matrices A and C.'
     WRITE( NOUT, FMT = 9999 ) &
                 'N       : The number of columns in the '// &
                 'matrices B and C.'
     WRITE( NOUT, FMT = 9999 ) &
                 'K       : The number of rows of B '// &
                 'and the number of columns of A.'
     WRITE( NOUT, FMT = 9999 ) &
                 'NB      : The size of the square blocks the'// &
                 ' matrices A, B and C are split into.'
     WRITE( NOUT, FMT = 9999 ) &
                 'P       : The number of process rows.'
     WRITE( NOUT, FMT = 9999 ) &
                 'Q       : The number of process columns.'
     WRITE( NOUT, FMT = * )
     WRITE( NOUT, FMT = 9999 ) &
                 'The following parameter values will be used:'
     WRITE( NOUT, FMT = 9998 ) 'M    ', M
     WRITE( NOUT, FMT = 9998 ) 'N    ', N
     WRITE( NOUT, FMT = 9998 ) 'K    ', K
     WRITE( NOUT, FMT = 9998 ) 'NB   ', NB
     WRITE( NOUT, FMT = 9998 ) 'P    ', NPROW
     WRITE( NOUT, FMT = 9998 ) 'Q    ', NPCOL
     WRITE( NOUT, FMT = * )
!
  else
!
!        If underlying system needs additional set up, do it now
!
     if ( NPROCS < 1 ) &
        CALL BLACS_SETUP( IAM, NPROCS )
!
!        Temporarily define blacs grid to include all processes so
!        information can be broadcast to all processes
!
     CALL BLACS_GET( -1, 0, ICTXT )
     CALL BLACS_GRIDINIT( ICTXT, 'Row-major', 1, NPROCS )

     CALL IGEBR2D( ICTXT, 'All', ' ', 6, 1, WORK, 6, 0, 0 )
     M     = WORK( 1 )
     N     = WORK( 2 )
     K     = WORK( 3 )
     NB    = WORK( 4 )
     NPROW = WORK( 5 )
     NPCOL = WORK( 6 )

  end if

  CALL BLACS_GRIDEXIT( ICTXT )

  RETURN

   20 WRITE( NOUT, FMT = 9997 )
  CLOSE( NIN )
  if ( NOUT /= 6 .AND. NOUT /= 0 ) &
     CLOSE( NOUT )
  CALL BLACS_ABORT( ICTXT, 1 )

  STOP

 9999 FORMAT( A )
 9998 FORMAT( 2X, A5, '   :        ', I6 )
 9997 FORMAT( ' Illegal input in file ',40A,'.  Aborting run.' )
end
SUBROUTINE LADD( J, K, I )
!
!*****************************************************************************80
!
!  -- ScaLAPACK routine (version 1.0) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     February 28, 1995
!
!     .. Array Arguments ..
  INTEGER            I(2), J(2), K(2)
!     ..
!
  INTEGER            IPOW16, IPOW15
  PARAMETER        ( IPOW16=2**16, IPOW15=2**15 )
!
  I(1) = MOD( K(1)+J(1), IPOW16 )
  I(2) = MOD( (K(1)+J(1)) / IPOW16+K(2)+J(2), IPOW15 )

  return
END
SUBROUTINE LMUL( K, J, I )
!
!*****************************************************************************80
!
!  -- ScaLAPACK routine (version 1.0) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     February 28, 1995
!
!     .. Array Arguments ..
  INTEGER            I(2), J(2), K(2)
!     ..
!     .. Parameters ..
  INTEGER            IPOW15, IPOW16, IPOW30
  PARAMETER        ( IPOW15=2**15, IPOW16=2**16, IPOW30=2**30 )
!     ..
!     .. Local Scalars ..
  INTEGER            KT, LT
!
  KT   = K(1)*J(1)
  if ( KT.LT.0 ) KT = (KT+IPOW30) + IPOW30
  I(1) = MOD(KT,IPOW16)
  LT   = K(1)*J(2) + K(2)*J(1)
  if ( LT.LT.0 ) LT = (LT+IPOW30) + IPOW30
  KT   = KT/IPOW16 + LT
  if ( KT.LT.0 ) KT = (KT+IPOW30) + IPOW30
  I(2) = MOD( KT, IPOW15 )

  RETURN
END
SUBROUTINE XJUMPM( JUMPM, MULT, IADD, IRANN, IRANM, IAM, ICM )
!
!*****************************************************************************80
!
!  -- ScaLAPACK routine (version 1.0) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     February 28, 1995
!
!     .. Scalar Arguments ..
  INTEGER            JUMPM
!     ..
!     .. Array Arguments ..
  INTEGER            IADD(2), IAM(2), ICM(2), IRANM(2), IRANN(2)
  INTEGER            MULT(2)
!     ..
!
  INTEGER            I
!     ..
!     .. Local Arrays ..
  INTEGER            J(2)
!     ..
!     .. External Subroutines ..
  EXTERNAL           LADD, LMUL
!
  if ( JUMPM > 0 ) THEN

     DO I = 1, 2
        IAM(I) = MULT(I)
        ICM(I) = IADD(I)
     end do

     DO I = 1, JUMPM-1
        CALL LMUL( IAM, MULT, J )
        IAM(1) = J(1)
        IAM(2) = J(2)
        CALL LMUL( ICM, MULT, J )
        CALL LADD( IADD, J, ICM )
     end do

     CALL LMUL( IRANN, IAM, J )
     CALL LADD( J, ICM, IRANM )

  else

     IRANM(1) = IRANN(1)
     IRANM(2) = IRANN(2)

  end if

  return
END
SUBROUTINE SETRAN( IRAN, IA, IC )
!
!*****************************************************************************80
!
!  -- ScaLAPACK routine (version 1.0) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     February 28, 1995
!
!     .. Array Arguments ..
  INTEGER            IA(2),  IC(2), IRAN(2)

!     .. Local Scalars ..
  INTEGER            I
!     ..
!     .. Local Arrays ..
  INTEGER            IAS(2),  ICS(2), IRAND(2)
!     ..
!     .. Common Blocks ..
  COMMON /RANCOM/    IRAND, IAS, ICS
  SAVE   /RANCOM/
!     ..
!     .. Executable Statements ..
!
  DO I = 1, 2
     IRAND(I) = IRAN(I)
     IAS(I)   = IA(I)
     ICS(I)   = IC(I)
  end do

  RETURN
END
SUBROUTINE JUMPIT( MULT, IADD, IRANN, IRANM )
!
!*****************************************************************************80
!
!  -- ScaLAPACK routine (version 1.0) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     February 28, 1995
!
!     .. Array Arguments ..
  INTEGER            IADD(2), IRANM(2), IRANN(2), MULT(2)
!     ..
  INTEGER            IAS(2), ICS(2), IRAND(2), J(2)
!     ..
!     .. External Subroutines ..
  EXTERNAL           LADD, LMUL
!     ..
!     .. Common Blocks ..
  COMMON /RANCOM/    IRAND, IAS, ICS
  SAVE   /RANCOM/
!     ..
!     .. Executable Statements ..
!
  CALL LMUL( IRANN, MULT, J )
  CALL LADD( J, IADD, IRANM )

  IRAND(1) = IRANM(1)
  IRAND(2) = IRANM(2)

  return
end
REAL FUNCTION PSRAND( IDUMM )
!
!*****************************************************************************80
!
!  -- ScaLAPACK routine (version 1.0) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     February 28, 1995
!
!     .. Scalar Arguments ..
  INTEGER            IDUMM
!     ..
!
  REAL               DIVFAC, POW16
  PARAMETER          ( DIVFAC=2.147483648E+9, POW16=6.5536E+4 )
!     ..
!     .. Local Arrays ..
  INTEGER            J( 2 )
!     ..
!     .. External Subroutines ..
  EXTERNAL           LADD, LMUL
!     ..
!     .. Intrinsic Functions ..
  INTRINSIC          REAL
!     ..
!     .. Common Blocks ..
  INTEGER            IAS(2), ICS(2), IRAND(2)
  COMMON /RANCOM/    IRAND, IAS, ICS
  SAVE   /RANCOM/
!     ..
!     .. Executable Statements ..
!
  PSRAND = ( REAL(IRAND(1)) + POW16 * REAL(IRAND(2)) ) / DIVFAC

  CALL LMUL( IRAND, IAS, J )
  CALL LADD( J, ICS, IRAND )

  return
end
DOUBLE PRECISION FUNCTION PDRAND( IDUMM )
!
!*****************************************************************************80
!
!  -- ScaLAPACK routine (version 1.0) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     February 28, 1995
!
!     .. Scalar Arguments ..
  INTEGER            IDUMM
!     ..
!     .. Parameters ..
  DOUBLE PRECISION   DIVFAC, POW16
  PARAMETER          ( DIVFAC=2.147483648D+9, POW16=6.5536D+4 )
  INTEGER            J(2)
!     ..
!     .. External Subroutines ..
  EXTERNAL           LADD, LMUL
!     ..
!     .. Common Blocks ..
  INTEGER            IAS(2), ICS(2), IRAND(2)
  COMMON /RANCOM/    IRAND, IAS, ICS
  SAVE   /RANCOM/
!
  PDRAND = ( DBLE(IRAND(1)) + POW16 * DBLE(IRAND(2)) ) / DIVFAC

  CALL LMUL ( IRAND, IAS, J )

  CALL LADD ( J, ICS, IRAND )

  return
end
