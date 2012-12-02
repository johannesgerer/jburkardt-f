      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, &
                       N4 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
!     ..
!
!  Purpose
!  =======
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
!  =========
!
!  ISPEC   (input) INTEGER
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
!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper case or
!          lower case.
!
!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.
!
!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.
!
! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
!
!  Further Details
!  ===============
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
!      IF( NB<=1 ) NB = MAX( 1, N )
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..
!
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
!
!     Invalid value for ISPEC
!
      ILAENV = -1
      RETURN
!
  100 CONTINUE
!
!     Convert NAME to upper case if the first character is lower case.
!
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ==90 .OR. IZ==122 ) THEN
!
!        ASCII character set
!
         IF( IC>=97 .AND. IC<=122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC>=97 .AND. IC<=122 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
!
      ELSE IF( IZ==233 .OR. IZ==169 ) THEN
!
!        EBCDIC character set
!
         IF( ( IC>=129 .AND. IC<=137 ) .OR. &
             ( IC>=145 .AND. IC<=153 ) .OR. &
             ( IC>=162 .AND. IC<=169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC>=129 .AND. IC<=137 ) .OR. &
                   ( IC>=145 .AND. IC<=153 ) .OR. &
                   ( IC>=162 .AND. IC<=169 ) ) &
                  SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
!
      ELSE IF( IZ==218 .OR. IZ==250 ) THEN
!
!        Prime machines:  ASCII+128
!
         IF( IC>=225 .AND. IC<=250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC>=225 .AND. IC<=250 ) &
                  SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
!
      C1 = SUBNAM( 1:1 )
      SNAME = C1=='S' .OR. C1=='D'
      CNAME = C1=='C' .OR. C1=='Z'
      IF( .NOT.( CNAME .OR. SNAME ) ) &
         RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
!
      GO TO ( 110, 200, 300 ) ISPEC
!
  110 CONTINUE
!
!     ISPEC = 1:  block size
!
!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.
!
      NB = 1
!
      IF( C2=='GE' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR. &
                  C3=='QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3=='HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3=='BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3=='TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2=='PO' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2=='SY' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3=='TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3=='GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2=='HE' ) THEN
         IF( C3=='TRF' ) THEN
            NB = 64
         ELSE IF( C3=='TRD' ) THEN
            NB = 1
         ELSE IF( C3=='GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2=='OR' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2=='UN' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2=='GB' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2=='PB' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2<=64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2=='TR' ) THEN
         IF( C3=='TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2=='LA' ) THEN
         IF( C3=='UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2=='ST' ) THEN
         IF( C3=='EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
!
  200 CONTINUE
!
!     ISPEC = 2:  minimum block size
!
      NBMIN = 2
      IF( C2=='GE' ) THEN
         IF( C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR. &
             C3=='QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3=='HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3=='BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3=='TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2=='SY' ) THEN
         IF( C3=='TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF( SNAME .AND. C3=='TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2=='HE' ) THEN
         IF( C3=='TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2=='OR' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2=='UN' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 )=='M' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
!
  300 CONTINUE
!
!     ISPEC = 3:  crossover point
!
      NX = 0
      IF( C2=='GE' ) THEN
         IF( C3=='QRF' .OR. C3=='RQF' .OR. C3=='LQF' .OR. &
             C3=='QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3=='HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3=='BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2=='SY' ) THEN
         IF( SNAME .AND. C3=='TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2=='HE' ) THEN
         IF( C3=='TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2=='OR' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2=='UN' ) THEN
         IF( C3( 1:1 )=='G' ) THEN
            IF( C4=='QR' .OR. C4=='RQ' .OR. C4=='LQ' .OR. &
                C4=='QL' .OR. C4=='HR' .OR. C4=='TR' .OR. &
                C4=='BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
!
  400 CONTINUE
!
!     ISPEC = 4:  number of shifts (used by xHSEQR)
!
      ILAENV = 6
      RETURN
!
  500 CONTINUE
!
!     ISPEC = 5:  minimum column dimension (not used)
!
      ILAENV = 2
      RETURN
!
  600 CONTINUE
!
!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
!
  700 CONTINUE
!
!     ISPEC = 7:  number of processors (not used)
!
      ILAENV = 1
      RETURN
!
  800 CONTINUE
!
!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      ILAENV = 50
      RETURN
!
!     End of ILAENV
!
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          CA, CB
!     ..
!
!  Purpose
!  =======
!
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of
!  case.
!
!  Arguments
!  =========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!          CA and CB specify the single characters to be compared.
!
! =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
!     ..
!     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
!     ..
!     .. Executable Statements ..
!
!     Test if the characters are equal
!
      LSAME = CA==CB
      IF( LSAME ) &
         RETURN
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
      IF( ZCODE==90 .OR. ZCODE==122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA>=97 .AND. INTA<=122 ) INTA = INTA - 32
         IF( INTB>=97 .AND. INTB<=122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE==233 .OR. ZCODE==169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA>=129 .AND. INTA<=137 .OR. &
             INTA>=145 .AND. INTA<=153 .OR. &
             INTA>=162 .AND. INTA<=169 ) INTA = INTA + 64
         IF( INTB>=129 .AND. INTB<=137 .OR. &
             INTB>=145 .AND. INTB<=153 .OR. &
             INTB>=162 .AND. INTB<=169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE==218 .OR. ZCODE==250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA>=225 .AND. INTA<=250 ) INTA = INTA - 32
         IF( INTB>=225 .AND. INTB<=250 ) INTB = INTB - 32
      END IF
      LSAME = INTA==INTB
!
!     RETURN
!
!     End of LSAME
!
      END
      LOGICAL          FUNCTION LSAMEN( N, CA, CB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*( * )    CA, CB
      INTEGER            N
!     ..
!
!  Purpose
!  =======
!
!  LSAMEN  tests if the first N letters of CA are the same as the
!  first N letters of CB, regardless of case.
!  LSAMEN returns .TRUE. if CA and CB are equivalent except for case
!  and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )
!  or LEN( CB ) is less than N.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of characters in CA and CB to be compared.
!
!  CA      (input) CHARACTER*(*)
!  CB      (input) CHARACTER*(*)
!          CA and CB specify two character strings of length at least N.
!          Only the first N characters of each string will be accessed.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          LEN
!     ..
!     .. Executable Statements ..
!
      LSAMEN = .FALSE.
      IF( LEN( CA )<N .OR. LEN( CB )<N ) &
         GO TO 20
!
!     Do for each character in the two strings.
!
      DO 10 I = 1, N
!
!        Test if the characters are equal using LSAME.
!
         IF( .NOT.LSAME( CA( I: I ), CB( I: I ) ) ) &
            GO TO 20
!
   10 CONTINUE
      LSAMEN = .TRUE.
!
   20 CONTINUE
      RETURN
!
!     End of LSAMEN
!
      END
      REAL             FUNCTION SCSUM1( N, CX, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
!     ..
!     .. Array Arguments ..
      COMPLEX            CX( * )
!     ..
!
!  Purpose
!  =======
!
!  SCSUM1 takes the sum of the absolute values of a complex
!  vector and returns a single precision result.
!
!  Based on SCASUM from the Level 1 BLAS.
!  The change is to use the 'genuine' absolute value.
!
!  Contributed by Nick Higham for use with CLACON.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements in the vector CX.
!
!  CX      (input) COMPLEX array, dimension (N)
!          The vector whose elements will be summed.
!
!  INCX    (input) INTEGER
!          The spacing between successive values of CX.  INCX > 0.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, NINCX
      REAL               STEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      SCSUM1 = 0.0E0
      STEMP = 0.0E0
      IF( N<=0 ) &
         RETURN
      IF( INCX==1 ) &
         GO TO 20
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
   20 CONTINUE
      DO 30 I = 1, N
!
!        NEXT LINE MODIFIED.
!
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      SCSUM1 = STEMP
      RETURN
!
!     End of SCSUM1
!
      END
      SUBROUTINE SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  SGBTF2 computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the unblocked version of the algorithm, calling Level 2 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) REAL array, dimension (LDAB,N)
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
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, JP, JU, KM, KV
!     ..
!     .. External Functions ..
      INTEGER            ISAMAX
      EXTERNAL           ISAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in.
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M<0 ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      ELSE IF( KL<0 ) THEN
         INFO = -3
      ELSE IF( KU<0 ) THEN
         INFO = -4
      ELSE IF( LDAB<KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SGBTF2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M==0 .OR. N==0 ) &
         RETURN
!
!     Gaussian elimination with partial pivoting
!
!     Set fill-in elements in columns KU+2 to KV to zero.
!
      DO 20 J = KU + 2, MIN( KV, N )
         DO 10 I = KV - J + 2, KL
            AB( I, J ) = ZERO
   10    CONTINUE
   20 CONTINUE
!
!     JU is the index of the last column affected by the current stage
!     of the factorization.
!
      JU = 1
!
      DO 40 J = 1, MIN( M, N )
!
!        Set fill-in elements in column J+KV to zero.
!
         IF( J+KV<=N ) THEN
            DO 30 I = 1, KL
               AB( I, J+KV ) = ZERO
   30       CONTINUE
         END IF
!
!        Find pivot and test for singularity. KM is the number of
!        subdiagonal elements in the current column.
!
         KM = MIN( KL, M-J )
         JP = ISAMAX( KM+1, AB( KV+1, J ), 1 )
         IPIV( J ) = JP + J - 1
         IF( AB( KV+JP, J )/=ZERO ) THEN
            JU = MAX( JU, MIN( J+KU+JP-1, N ) )
!
!           Apply interchange to columns J to JU.
!
            IF( JP/=1 ) &
               CALL SSWAP( JU-J+1, AB( KV+JP, J ), LDAB-1, &
                           AB( KV+1, J ), LDAB-1 )
!
            IF( KM>0 ) THEN
!
!              Compute multipliers.
!
               CALL SSCAL( KM, ONE / AB( KV+1, J ), AB( KV+2, J ), 1 )
!
!              Update trailing submatrix within the band.
!
               IF( JU>J ) &
                  CALL SGER( KM, JU-J, -ONE, AB( KV+2, J ), 1, &
                             AB( KV, J+1 ), LDAB-1, AB( KV+1, J+1 ), &
                             LDAB-1 )
            END IF
         ELSE
!
!           If pivot is zero, set INFO to the index of the pivot
!           unless a zero pivot has already been found.
!
            IF( INFO==0 ) &
               INFO = J
         END IF
   40 CONTINUE
      RETURN
!
!     End of SGBTF2
!
      END
      SUBROUTINE SGBTRF( M, N, KL, KU, AB, LDAB, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * )
!     ..
!
!  Purpose
!  =======
!
!  SGBTRF computes an LU factorization of a real m-by-n band matrix A
!  using partial pivoting with row interchanges.
!
!  This is the blocked version of the algorithm, calling Level 3 BLAS.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  AB      (input/output) REAL array, dimension (LDAB,N)
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
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.
!
!  Further Details
!  ===============
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NBMAX, LDWORK
      PARAMETER          ( NBMAX = 64, LDWORK = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I2, I3, II, IP, J, J2, J3, JB, JJ, JM, JP, &
                         JU, K2, KM, KV, NB, NW
      REAL               TEMP
!     ..
!     .. Local Arrays ..
      REAL               WORK13( LDWORK, NBMAX ), &
                         WORK31( LDWORK, NBMAX )
!     ..
!     .. External Functions ..
      INTEGER            ILAENV, ISAMAX
      EXTERNAL           ILAENV, ISAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SCOPY, SGBTF2, SGEMM, SGER, SLASWP, SSCAL, &
                         SSWAP, STRSM, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     KV is the number of superdiagonals in the factor U, allowing for
!     fill-in
!
      KV = KU + KL
!
!     Test the input parameters.
!
      INFO = 0
      IF( M<0 ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      ELSE IF( KL<0 ) THEN
         INFO = -3
      ELSE IF( KU<0 ) THEN
         INFO = -4
      ELSE IF( LDAB<KL+KV+1 ) THEN
         INFO = -6
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SGBTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M==0 .OR. N==0 ) &
         RETURN
!
!     Determine the block size for this environment
!
      NB = ILAENV( 1, 'SGBTRF', ' ', M, N, KL, KU )
!
!     The block size must not exceed the limit set by the size of the
!     local arrays WORK13 and WORK31.
!
      NB = MIN( NB, NBMAX )
!
      IF( NB<=1 .OR. NB>KL ) THEN
!
!        Use unblocked code
!
         CALL SGBTF2( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      ELSE
!
!        Use blocked code
!
!        Zero the superdiagonal elements of the work array WORK13
!
         DO 20 J = 1, NB
            DO 10 I = 1, J - 1
               WORK13( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
!
!        Zero the subdiagonal elements of the work array WORK31
!
         DO 40 J = 1, NB
            DO 30 I = J + 1, NB
               WORK31( I, J ) = ZERO
   30       CONTINUE
   40    CONTINUE
!
!        Gaussian elimination with partial pivoting
!
!        Set fill-in elements in columns KU+2 to KV to zero
!
         DO 60 J = KU + 2, MIN( KV, N )
            DO 50 I = KV - J + 2, KL
               AB( I, J ) = ZERO
   50       CONTINUE
   60    CONTINUE
!
!        JU is the index of the last column affected by the current
!        stage of the factorization
!
         JU = 1
!
         DO 180 J = 1, MIN( M, N ), NB
            JB = MIN( NB, MIN( M, N )-J+1 )
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
            I2 = MIN( KL-JB, M-J-JB+1 )
            I3 = MIN( JB, M-J-KL+1 )
!
!           J2 and J3 are computed after JU has been updated.
!
!           Factorize the current block of JB columns
!
            DO 80 JJ = J, J + JB - 1
!
!              Set fill-in elements in column JJ+KV to zero
!
               IF( JJ+KV<=N ) THEN
                  DO 70 I = 1, KL
                     AB( I, JJ+KV ) = ZERO
   70             CONTINUE
               END IF
!
!              Find pivot and test for singularity. KM is the number of
!              subdiagonal elements in the current column.
!
               KM = MIN( KL, M-JJ )
               JP = ISAMAX( KM+1, AB( KV+1, JJ ), 1 )
               IPIV( JJ ) = JP + JJ - J
               IF( AB( KV+JP, JJ )/=ZERO ) THEN
                  JU = MAX( JU, MIN( JJ+KU+JP-1, N ) )
                  IF( JP/=1 ) THEN
!
!                    Apply interchange to columns J to J+JB-1
!
                     IF( JP+JJ-1<J+KL ) THEN
!
                        CALL SSWAP( JB, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    AB( KV+JP+JJ-J, J ), LDAB-1 )
                     ELSE
!
!                       The interchange affects columns J to JJ-1 of A31
!                       which are stored in the work array WORK31
!
                        CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                    WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                        CALL SSWAP( J+JB-JJ, AB( KV+1, JJ ), LDAB-1, &
                                    AB( KV+JP, JJ ), LDAB-1 )
                     END IF
                  END IF
!
!                 Compute multipliers
!
                  CALL SSCAL( KM, ONE / AB( KV+1, JJ ), AB( KV+2, JJ ), &
                              1 )
!
!                 Update trailing submatrix within the band and within
!                 the current block. JM is the index of the last column
!                 which needs to be updated.
!
                  JM = MIN( JU, J+JB-1 )
                  IF( JM>JJ ) &
                     CALL SGER( KM, JM-JJ, -ONE, AB( KV+2, JJ ), 1, &
                                AB( KV, JJ+1 ), LDAB-1, &
                                AB( KV+1, JJ+1 ), LDAB-1 )
               ELSE
!
!                 If pivot is zero, set INFO to the index of the pivot
!                 unless a zero pivot has already been found.
!
                  IF( INFO==0 ) &
                     INFO = JJ
               END IF
!
!              Copy current column of A31 into the work array WORK31
!
               NW = MIN( JJ-J+1, I3 )
               IF( NW>0 ) &
                  CALL SCOPY( NW, AB( KV+KL+1-JJ+J, JJ ), 1, &
                              WORK31( 1, JJ-J+1 ), 1 )
   80       CONTINUE
            IF( J+JB<=N ) THEN
!
!              Apply the row interchanges to the other blocks.
!
               J2 = MIN( JU-J+1, KV ) - JB
               J3 = MAX( 0, JU-J-KV+1 )
!
!              Use SLASWP to apply the row interchanges to A12, A22, and
!              A32.
!
               CALL SLASWP( J2, AB( KV+1-JB, J+JB ), LDAB-1, 1, JB, &
                            IPIV( J ), 1 )
!
!              Adjust the pivot indices.
!
               DO 90 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
   90          CONTINUE
!
!              Apply the row interchanges to A13, A23, and A33
!              columnwise.
!
               K2 = J - 1 + JB + J2
               DO 110 I = 1, J3
                  JJ = K2 + I
                  DO 100 II = J + I - 1, J + JB - 1
                     IP = IPIV( II )
                     IF( IP/=II ) THEN
                        TEMP = AB( KV+1+II-JJ, JJ )
                        AB( KV+1+II-JJ, JJ ) = AB( KV+1+IP-JJ, JJ )
                        AB( KV+1+IP-JJ, JJ ) = TEMP
                     END IF
  100             CONTINUE
  110          CONTINUE
!
!              Update the relevant part of the trailing submatrix
!
               IF( J2>0 ) THEN
!
!                 Update A12
!
                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J2, ONE, AB( KV+1, J ), LDAB-1, &
                              AB( KV+1-JB, J+JB ), LDAB-1 )
!
                  IF( I2>0 ) THEN
!
!                    Update A22
!
                     CALL SGEMM( 'No transpose', 'No transpose', I2, J2, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+1, J+JB ), LDAB-1 )
                  END IF
!
                  IF( I3>0 ) THEN
!
!                    Update A32
!
                     CALL SGEMM( 'No transpose', 'No transpose', I3, J2, &
                                 JB, -ONE, WORK31, LDWORK, &
                                 AB( KV+1-JB, J+JB ), LDAB-1, ONE, &
                                 AB( KV+KL+1-JB, J+JB ), LDAB-1 )
                  END IF
               END IF
!
               IF( J3>0 ) THEN
!
!                 Copy the lower triangle of A13 into the work array
!                 WORK13
!
                  DO 130 JJ = 1, J3
                     DO 120 II = JJ, JB
                        WORK13( II, JJ ) = AB( II-JJ+1, JJ+J+KV-1 )
  120                CONTINUE
  130             CONTINUE
!
!                 Update A13 in the work array
!
                  CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', &
                              JB, J3, ONE, AB( KV+1, J ), LDAB-1, &
                              WORK13, LDWORK )
!
                  IF( I2>0 ) THEN
!
!                    Update A23
!
                     CALL SGEMM( 'No transpose', 'No transpose', I2, J3, &
                                 JB, -ONE, AB( KV+1+JB, J ), LDAB-1, &
                                 WORK13, LDWORK, ONE, AB( 1+JB, J+KV ), &
                                 LDAB-1 )
                  END IF
!
                  IF( I3>0 ) THEN
!
!                    Update A33
!
                     CALL SGEMM( 'No transpose', 'No transpose', I3, J3, &
                                 JB, -ONE, WORK31, LDWORK, WORK13, &
                                 LDWORK, ONE, AB( 1+KL, J+KV ), LDAB-1 )
                  END IF
!
!                 Copy the lower triangle of A13 back into place
!
                  DO 150 JJ = 1, J3
                     DO 140 II = JJ, JB
                        AB( II-JJ+1, JJ+J+KV-1 ) = WORK13( II, JJ )
  140                CONTINUE
  150             CONTINUE
               END IF
            ELSE
!
!              Adjust the pivot indices.
!
               DO 160 I = J, J + JB - 1
                  IPIV( I ) = IPIV( I ) + J - 1
  160          CONTINUE
            END IF
!
!           Partially undo the interchanges in the current block to
!           restore the upper triangular form of A31 and copy the upper
!           triangle of A31 back into place
!
            DO 170 JJ = J + JB - 1, J, -1
               JP = IPIV( JJ ) - JJ + 1
               IF( JP/=1 ) THEN
!
!                 Apply interchange to columns J to JJ-1
!
                  IF( JP+JJ-1<J+KL ) THEN
!
!                    The interchange does not affect A31
!
                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 AB( KV+JP+JJ-J, J ), LDAB-1 )
                  ELSE
!
!                    The interchange does affect A31
!
                     CALL SSWAP( JJ-J, AB( KV+1+JJ-J, J ), LDAB-1, &
                                 WORK31( JP+JJ-J-KL, 1 ), LDWORK )
                  END IF
               END IF
!
!              Copy the current column of A31 back into place
!
               NW = MIN( I3, JJ-J+1 )
               IF( NW>0 ) &
                  CALL SCOPY( NW, WORK31( 1, JJ-J+1 ), 1, &
                              AB( KV+KL+1-JJ+J, JJ ), 1 )
  170       CONTINUE
  180    CONTINUE
      END IF
!
      RETURN
!
!     End of SGBTRF
!
      END
      SUBROUTINE SGBTRS( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               AB( LDAB, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  SGBTRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general band matrix A using the LU factorization computed
!  by SGBTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  KL      (input) INTEGER
!          The number of subdiagonals within the band of A.  KL >= 0.
!
!  KU      (input) INTEGER
!          The number of superdiagonals within the band of A.  KU >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AB      (input) REAL array, dimension (LDAB,N)
!          Details of the LU factorization of the band matrix A, as
!          computed by SGBTRF.  U is stored as an upper triangular band
!          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!          the multipliers used during the factorization are stored in
!          rows KL+KU+2 to 2*KL+KU+1.
!
!  LDAB    (input) INTEGER
!          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= N, row i of the matrix was
!          interchanged with row IPIV(i).
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            I, J, KD, L, LM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER, SSWAP, STBSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
          LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      ELSE IF( KL<0 ) THEN
         INFO = -3
      ELSE IF( KU<0 ) THEN
         INFO = -4
      ELSE IF( NRHS<0 ) THEN
         INFO = -5
      ELSE IF( LDAB<( 2*KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB<MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SGBTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N==0 .OR. NRHS==0 ) &
         RETURN
!
      KD = KU + KL + 1
      LNOTI = KL>0
!
      IF( NOTRAN ) THEN
!
!        Solve  A*X = B.
!
!        Solve L*X = B, overwriting B with X.
!
!        L is represented as a product of permutations and unit lower
!        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
!        where each transformation L(i) is a rank-one modification of
!        the identity matrix.
!
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L/=J ) &
                  CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
               CALL SGER( LM, NRHS, -ONE, AB( KD+1, J ), 1, B( J, 1 ), &
                          LDB, B( J+1, 1 ), LDB )
   10       CONTINUE
         END IF
!
         DO 20 I = 1, NRHS
!
!           Solve U*X = B, overwriting B with X.
!
            CALL STBSV( 'Upper', 'No transpose', 'Non-unit', N, KL+KU, &
                        AB, LDAB, B( 1, I ), 1 )
   20    CONTINUE
!
      ELSE
!
!        Solve A'*X = B.
!
         DO 30 I = 1, NRHS
!
!           Solve U'*X = B, overwriting B with X.
!
            CALL STBSV( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB, &
                        LDAB, B( 1, I ), 1 )
   30    CONTINUE
!
!        Solve L'*X = B, overwriting B with X.
!
         IF( LNOTI ) THEN
            DO 40 J = N - 1, 1, -1
               LM = MIN( KL, N-J )
               CALL SGEMV( 'Transpose', LM, NRHS, -ONE, B( J+1, 1 ), &
                           LDB, AB( KD+1, J ), 1, ONE, B( J, 1 ), LDB )
               L = IPIV( J )
               IF( L/=J ) &
                  CALL SSWAP( NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB )
   40       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of SGBTRS
!
      END
      SUBROUTINE SGEQR2( M, N, A, LDA, TAU, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SGEQR2 computes a QR factorization of a real m by n matrix A:
!  A = Q * R.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the m by n matrix A.
!          On exit, the elements on and above the diagonal of the array
!          contain the min(m,n) by n upper trapezoidal matrix R (R is
!          upper triangular if m >= n); the elements below the diagonal,
!          with the array TAU, represent the orthogonal matrix Q as a
!          product of elementary reflectors (see Further Details).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  TAU     (output) REAL array, dimension (min(M,N))
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace) REAL array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
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
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, K
      REAL               AII
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, SLARFG, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      IF( M<0 ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      ELSE IF( LDA<MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SGEQR2', -INFO )
         RETURN
      END IF
!
      K = MIN( M, N )
!
      DO 10 I = 1, K
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
         CALL SLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1, &
                      TAU( I ) )
         IF( I<N ) THEN
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
            AII = A( I, I )
            A( I, I ) = ONE
            CALL SLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ), &
                        A( I, I+1 ), LDA, WORK )
            A( I, I ) = AII
         END IF
   10 CONTINUE
      RETURN
!
!     End of SGEQR2
!
      END
      SUBROUTINE SGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               D( * ), DL( * ), DU( * ), DU2( * )
!     ..
!
!  Purpose
!  =======
!
!  SGTTRF computes an LU factorization of a real tridiagonal matrix A
!  using elimination with partial pivoting and row interchanges.
!
!  The factorization has the form
!     A = L * U
!  where L is a product of permutation and unit lower bidiagonal
!  matrices and U is upper triangular with nonzeros in only the main
!  diagonal and first two superdiagonals.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  DL      (input/output) REAL array, dimension (N-1)
!          On entry, DL must contain the (n-1) subdiagonal elements of
!          A.
!          On exit, DL is overwritten by the (n-1) multipliers that
!          define the matrix L from the LU factorization of A.
!
!  D       (input/output) REAL array, dimension (N)
!          On entry, D must contain the diagonal elements of A.
!          On exit, D is overwritten by the n diagonal elements of the
!          upper triangular matrix U from the LU factorization of A.
!
!  DU      (input/output) REAL array, dimension (N-1)
!          On entry, DU must contain the (n-1) superdiagonal elements
!          of A.
!          On exit, DU is overwritten by the (n-1) elements of the first
!          superdiagonal of U.
!
!  DU2     (output) REAL array, dimension (N-2)
!          On exit, DU2 is overwritten by the (n-2) elements of the
!          second superdiagonal of U.
!
!  IPIV    (output) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= n, row i of the matrix was
!          interchanged with row IPIV(i).  IPIV(i) will always be either
!          i or i+1; IPIV(i) = i indicates a row interchange was not
!          required.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      REAL               FACT, TEMP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      IF( N<0 ) THEN
         INFO = -1
         CALL XERBLA( 'SGTTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N==0 ) &
         RETURN
!
!     Initialize IPIV(i) = i
!
      DO 10 I = 1, N
         IPIV( I ) = I
   10 CONTINUE
!
      DO 20 I = 1, N - 1
         IF( DL( I )==ZERO ) THEN
!
!           Subdiagonal is zero, no elimination is required.
!
            IF( D( I )==ZERO .AND. INFO==0 ) &
               INFO = I
            IF( I<N-1 ) &
               DU2( I ) = ZERO
         ELSE IF( ABS( D( I ) )>=ABS( DL( I ) ) ) THEN
!
!           No row interchange required, eliminate DL(I)
!
            FACT = DL( I ) / D( I )
            DL( I ) = FACT
            D( I+1 ) = D( I+1 ) - FACT*DU( I )
            IF( I<N-1 ) &
               DU2( I ) = ZERO
         ELSE
!
!           Interchange rows I and I+1, eliminate DL(I)
!
            FACT = D( I ) / DL( I )
            D( I ) = DL( I )
            DL( I ) = FACT
            TEMP = DU( I )
            DU( I ) = D( I+1 )
            D( I+1 ) = TEMP - FACT*D( I+1 )
            IF( I<N-1 ) THEN
               DU2( I ) = DU( I+1 )
               DU( I+1 ) = -FACT*DU( I+1 )
            END IF
            IPIV( I ) = IPIV( I ) + 1
         END IF
   20 CONTINUE
      IF( D( N )==ZERO .AND. INFO==0 ) THEN
         INFO = N
         RETURN
      END IF
!
      RETURN
!
!     End of SGTTRF
!
      END
      SUBROUTINE SGTTRS( TRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
!     ..
!
!  Purpose
!  =======
!
!  SGTTRS solves one of the systems of equations
!     A*X = B  or  A'*X = B,
!  with a tridiagonal matrix A using the LU factorization computed
!  by SGTTRF.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER
!          Specifies the form of the system of equations:
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  DL      (input) REAL array, dimension (N-1)
!          The (n-1) multipliers that define the matrix L from the
!          LU factorization of A.
!
!  D       (input) REAL array, dimension (N)
!          The n diagonal elements of the upper triangular matrix U from
!          the LU factorization of A.
!
!  DU      (input) REAL array, dimension (N-1)
!          The (n-1) elements of the first superdiagonal of U.
!
!  DU2     (input) REAL array, dimension (N-2)
!          The (n-2) elements of the second superdiagonal of U.
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices; for 1 <= i <= n, row i of the matrix was
!          interchanged with row IPIV(i).  IPIV(i) will always be either
!          i or i+1; IPIV(i) = i indicates a row interchange was not
!          required.
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, B is overwritten by the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            NOTRAN
      INTEGER            I, J
      REAL               TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT. &
          LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      ELSE IF( NRHS<0 ) THEN
         INFO = -3
      ELSE IF( LDB<MAX( N, 1 ) ) THEN
         INFO = -10
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SGTTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N==0 .OR. NRHS==0 ) &
         RETURN
!
      IF( NOTRAN ) THEN
!
!        Solve A*X = B using the LU factorization of A,
!        overwriting each right hand side vector with its solution.
!
         DO 30 J = 1, NRHS
!
!           Solve L*x = b.
!
            DO 10 I = 1, N - 1
               IF( IPIV( I )==I ) THEN
                  B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
               ELSE
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - DL( I )*B( I, J )
               END IF
   10       CONTINUE
!
!           Solve U*x = b.
!
            B( N, J ) = B( N, J ) / D( N )
            IF( N>1 ) &
               B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / &
                             D( N-1 )
            DO 20 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )* &
                           B( I+2, J ) ) / D( I )
   20       CONTINUE
   30    CONTINUE
      ELSE
!
!        Solve A' * X = B.
!
         DO 60 J = 1, NRHS
!
!           Solve U'*x = b.
!
            B( 1, J ) = B( 1, J ) / D( 1 )
            IF( N>1 ) &
               B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
            DO 40 I = 3, N
               B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-DU2( I-2 )* &
                           B( I-2, J ) ) / D( I )
   40       CONTINUE
!
!           Solve L'*x = b.
!
            DO 50 I = N - 1, 1, -1
               IF( IPIV( I )==I ) THEN
                  B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
               ELSE
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                  B( I, J ) = TEMP
               END IF
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     End of SGTTRS
!
      END
      SUBROUTINE SLABAD( SMALL, LARGE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               LARGE, SMALL
!     ..
!
!  Purpose
!  =======
!
!  SLABAD takes as input the values computed by SLAMCH for underflow and
!  overflow, and returns the square root of each of these values if the
!  log of LARGE is sufficiently large.  This subroutine is intended to
!  identify machines with a large exponent range, such as the Crays, and
!  redefine the underflow and overflow limits to be the square roots of
!  the values computed by SLAMCH.  This subroutine is needed because
!  SLAMCH does not compensate for poor arithmetic in the upper half of
!  the exponent range, as is found on a Cray.
!
!  Arguments
!  =========
!
!  SMALL   (input/output) REAL
!          On entry, the underflow threshold as computed by SLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of SMALL, otherwise unchanged.
!
!  LARGE   (input/output) REAL
!          On entry, the overflow threshold as computed by SLAMCH.
!          On exit, if LOG10(LARGE) is sufficiently large, the square
!          root of LARGE, otherwise unchanged.
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF( LOG10( LARGE )>2000. ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
!
      RETURN
!
!     End of SLABAD
!
      END
      SUBROUTINE SLACON( N, V, X, ISGN, EST, KASE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            KASE, N
      REAL               EST
!     ..
!     .. Array Arguments ..
      INTEGER            ISGN( * )
      REAL               V( * ), X( * )
!     ..
!
!  Purpose
!  =======
!
!  SLACON estimates the 1-norm of a square, real matrix A.
!  Reverse communication is used for evaluating matrix-vector products.
!
!  Arguments
!  =========
!
!  N      (input) INTEGER
!         The order of the matrix.  N >= 1.
!
!  V      (workspace) REAL array, dimension (N)
!         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
!         (W is not returned).
!
!  X      (input/output) REAL array, dimension (N)
!         On an intermediate return, X should be overwritten by
!               A * X,   if KASE=1,
!               A' * X,  if KASE=2,
!         and SLACON must be re-called with all the other parameters
!         unchanged.
!
!  ISGN   (workspace) INTEGER array, dimension (N)
!
!  EST    (output) REAL
!         An estimate (a lower bound) for norm(A).
!
!  KASE   (input/output) INTEGER
!         On the initial call to SLACON, KASE should be 0.
!         On an intermediate return, KASE will be 1 or 2, indicating
!         whether X should be overwritten by A * X  or A' * X.
!         On the final return from SLACON, KASE will again be 0.
!
!  Further Details
!  ======= =======
!
!  Contributed by Nick Higham, University of Manchester.
!  Originally named SONEST, dated March 16, 1988.
!
!  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
!  a real or complex matrix, with applications to condition estimation",
!  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0, TWO = 2.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ITER, J, JLAST, JUMP
      REAL               ALTSGN, ESTOLD, TEMP
!     ..
!     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SASUM
      EXTERNAL           ISAMAX, SASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL           SCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, NINT, REAL, SIGN
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Executable Statements ..
!
      IF( KASE==0 ) THEN
         DO 10 I = 1, N
            X( I ) = ONE / REAL( N )
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
!
      GO TO ( 20, 40, 70, 110, 140 )JUMP
!
!     ................ ENTRY   (JUMP = 1)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
!
   20 CONTINUE
      IF( N==1 ) THEN
         V( 1 ) = X( 1 )
         EST = ABS( V( 1 ) )
!        ... QUIT
         GO TO 150
      END IF
      EST = SASUM( N, X, 1 )
!
      DO 30 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
!
!     ................ ENTRY   (JUMP = 2)
!     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
   40 CONTINUE
      J = ISAMAX( N, X, 1 )
      ITER = 2
!
!     MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
!
   50 CONTINUE
      DO 60 I = 1, N
         X( I ) = ZERO
   60 CONTINUE
      X( J ) = ONE
      KASE = 1
      JUMP = 3
      RETURN
!
!     ................ ENTRY   (JUMP = 3)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
   70 CONTINUE
      CALL SCOPY( N, X, 1, V, 1 )
      ESTOLD = EST
      EST = SASUM( N, V, 1 )
      DO 80 I = 1, N
         IF( NINT( SIGN( ONE, X( I ) ) )/=ISGN( I ) ) &
            GO TO 90
   80 CONTINUE
!     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
      GO TO 120
!
   90 CONTINUE
!     TEST FOR CYCLING.
      IF( EST<=ESTOLD ) &
         GO TO 120
!
      DO 100 I = 1, N
         X( I ) = SIGN( ONE, X( I ) )
         ISGN( I ) = NINT( X( I ) )
  100 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
!
!     ................ ENTRY   (JUMP = 4)
!     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X.
!
  110 CONTINUE
      JLAST = J
      J = ISAMAX( N, X, 1 )
      IF( ( X( JLAST )/=ABS( X( J ) ) ) .AND. ( ITER<ITMAX ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
!
!     ITERATION COMPLETE.  FINAL STAGE.
!
  120 CONTINUE
      ALTSGN = ONE
      DO 130 I = 1, N
         X( I ) = ALTSGN*( ONE+REAL( I-1 ) / REAL( N-1 ) )
         ALTSGN = -ALTSGN
  130 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
!
!     ................ ENTRY   (JUMP = 5)
!     X HAS BEEN OVERWRITTEN BY A*X.
!
  140 CONTINUE
      TEMP = TWO*( SASUM( N, X, 1 ) / REAL( 3*N ) )
      IF( TEMP>EST ) THEN
         CALL SCOPY( N, X, 1, V, 1 )
         EST = TEMP
      END IF
!
  150 CONTINUE
      KASE = 0
      RETURN
!
!     End of SLACON
!
      END
      SUBROUTINE SLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  SLACPY copies all or part of a two-dimensional matrix A to another
!  matrix B.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be copied to B.
!          = 'U':      Upper triangular part
!          = 'L':      Lower triangular part
!          Otherwise:  All of the matrix A
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!          or trapezoid is accessed; if UPLO = 'L', only the lower
!          triangle or trapezoid is accessed.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  B       (output) REAL array, dimension (LDB,N)
!          On exit, B = A in the locations specified by UPLO.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,M).
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
!
!     End of SLACPY
!
      END
      SUBROUTINE SLADIV( A, B, C, D, P, Q )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               A, B, C, D, P, Q
!     ..
!
!  Purpose
!  =======
!
!  SLADIV performs complex division in  real arithmetic
!
!                        a + i*b
!             p + i*q = ---------
!                        c + i*d
!
!  The algorithm is due to Robert L. Smith and can be found
!  in D. Knuth, The art of Computer Programming, Vol.2, p.195
!
!  Arguments
!  =========
!
!  A       (input) REAL
!  B       (input) REAL
!  C       (input) REAL
!  D       (input) REAL
!          The scalars a, b, c, and d in the above expression.
!
!  P       (output) REAL
!  Q       (output) REAL
!          The scalars p and q in the above expression.
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL               E, F
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( ABS( D )<ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
!
      RETURN
!
!     End of SLADIV
!
      END
      SUBROUTINE SLAE2( A, B, C, RT1, RT2 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               A, B, C, RT1, RT2
!     ..
!
!  Purpose
!  =======
!
!  SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, and RT2
!  is the eigenvalue of smaller absolute value.
!
!  Arguments
!  =========
!
!  A       (input) REAL
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) REAL
!          The (1,2) and (2,1) elements of the 2-by-2 matrix.
!
!  C       (input) REAL
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) REAL
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) REAL
!          The eigenvalue of smaller absolute value.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
!     ..
!     .. Local Scalars ..
      REAL               AB, ACMN, ACMX, ADF, DF, RT, SM, TB
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A )>ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF>AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF<AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
!
!        Includes case AB=ADF=0
!
         RT = AB*SQRT( TWO )
      END IF
      IF( SM<ZERO ) THEN
         RT1 = HALF*( SM-RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM>ZERO ) THEN
         RT1 = HALF*( SM+RT )
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         RT1 = HALF*RT
         RT2 = -HALF*RT
      END IF
      RETURN
!
!     End of SLAE2
!
      END
      SUBROUTINE SLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               A, B, C, CS1, RT1, RT2, SN1
!     ..
!
!  Purpose
!  =======
!
!  SLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition
!
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
!
!  Arguments
!  =========
!
!  A       (input) REAL
!          The (1,1) element of the 2-by-2 matrix.
!
!  B       (input) REAL
!          The (1,2) element and the conjugate of the (2,1) element of
!          the 2-by-2 matrix.
!
!  C       (input) REAL
!          The (2,2) element of the 2-by-2 matrix.
!
!  RT1     (output) REAL
!          The eigenvalue of larger absolute value.
!
!  RT2     (output) REAL
!          The eigenvalue of smaller absolute value.
!
!  CS1     (output) REAL
!  SN1     (output) REAL
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.
!
!  Further Details
!  ===============
!
!  RT1 is accurate to a few ulps barring over/underflow.
!
!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.
!
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.
!
!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               HALF
      PARAMETER          ( HALF = 0.5E0 )
!     ..
!     .. Local Scalars ..
      INTEGER            SGN1, SGN2
      REAL               AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, &
                         TB, TN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      SM = A + C
      DF = A - C
      ADF = ABS( DF )
      TB = B + B
      AB = ABS( TB )
      IF( ABS( A )>ABS( C ) ) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
      IF( ADF>AB ) THEN
         RT = ADF*SQRT( ONE+( AB / ADF )**2 )
      ELSE IF( ADF<AB ) THEN
         RT = AB*SQRT( ONE+( ADF / AB )**2 )
      ELSE
!
!        Includes case AB=ADF=0
!
         RT = AB*SQRT( TWO )
      END IF
      IF( SM<ZERO ) THEN
         RT1 = HALF*( SM-RT )
         SGN1 = -1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE IF( SM>ZERO ) THEN
         RT1 = HALF*( SM+RT )
         SGN1 = 1
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         RT2 = ( ACMX / RT1 )*ACMN - ( B / RT1 )*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         RT1 = HALF*RT
         RT2 = -HALF*RT
         SGN1 = 1
      END IF
!
!     Compute the eigenvector
!
      IF( DF>=ZERO ) THEN
         CS = DF + RT
         SGN2 = 1
      ELSE
         CS = DF - RT
         SGN2 = -1
      END IF
      ACS = ABS( CS )
      IF( ACS>AB ) THEN
         CT = -TB / CS
         SN1 = ONE / SQRT( ONE+CT*CT )
         CS1 = CT*SN1
      ELSE
         IF( AB==ZERO ) THEN
            CS1 = ONE
            SN1 = ZERO
         ELSE
            TN = -CS / TB
            CS1 = ONE / SQRT( ONE+TN*TN )
            SN1 = TN*CS1
         END IF
      END IF
      IF( SGN1==SGN2 ) THEN
         TN = CS1
         CS1 = -SN1
         SN1 = TN
      END IF
      RETURN
!
!     End of SLAEV2
!
      END
      SUBROUTINE SLAEXC( WANTQ, N, T, LDT, Q, LDQ, J1, N1, N2, WORK, &
                         INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            WANTQ
      INTEGER            INFO, J1, LDQ, LDT, N, N1, N2
!     ..
!     .. Array Arguments ..
      REAL               Q( LDQ, * ), T( LDT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLAEXC swaps adjacent diagonal blocks T11 and T22 of order 1 or 2 in
!  an upper quasi-triangular matrix T by an orthogonal similarity
!  transformation.
!
!  T must be in Schur canonical form, that is, block upper triangular
!  with 1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block
!  has its diagonal elemnts equal and its off-diagonal elements of
!  opposite sign.
!
!  Arguments
!  =========
!
!  WANTQ   (input) LOGICAL
!          = .TRUE. : accumulate the transformation in the matrix Q;
!          = .FALSE.: do not accumulate the transformation.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) REAL array, dimension (LDT,N)
!          On entry, the upper quasi-triangular matrix T, in Schur
!          canonical form.
!          On exit, the updated matrix T, again in Schur canonical form.
!
!  LDT     (input)  INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  Q       (input/output) REAL array, dimension (LDQ,N)
!          On entry, if WANTQ is .TRUE., the orthogonal matrix Q.
!          On exit, if WANTQ is .TRUE., the updated matrix Q.
!          If WANTQ is .FALSE., Q is not referenced.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q.
!          LDQ >= 1; and if WANTQ is .TRUE., LDQ >= N.
!
!  J1      (input) INTEGER
!          The index of the first row of the first block T11.
!
!  N1      (input) INTEGER
!          The order of the first block T11. N1 = 0, 1 or 2.
!
!  N2      (input) INTEGER
!          The order of the second block T22. N2 = 0, 1 or 2.
!
!  WORK    (workspace) REAL array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          = 1: the transformed matrix T would be too far from Schur
!               form; the blocks are not swapped and T and Q are
!               unchanged.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               TEN
      PARAMETER          ( TEN = 1.0E+1 )
      INTEGER            LDD, LDX
      PARAMETER          ( LDD = 4, LDX = 2 )
!     ..
!     .. Local Scalars ..
      INTEGER            IERR, J2, J3, J4, K, ND
      REAL               CS, DNORM, EPS, SCALE, SMLNUM, SN, T11, T22, &
                         T33, TAU, TAU1, TAU2, TEMP, THRESH, WI1, WI2, &
                         WR1, WR2, XNORM
!     ..
!     .. Local Arrays ..
      REAL               D( LDD, 4 ), U( 3 ), U1( 3 ), U2( 3 ), &
                         X( LDX, 2 )
!     ..
!     .. External Functions ..
      REAL               SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLACPY, SLANV2, SLARFG, SLARFX, SLARTG, SLASY2, &
                         SROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N==0 .OR. N1==0 .OR. N2==0 ) &
         RETURN
      IF( J1+N1>N ) &
         RETURN
!
      J2 = J1 + 1
      J3 = J1 + 2
      J4 = J1 + 3
!
      IF( N1==1 .AND. N2==1 ) THEN
!
!        Swap two 1-by-1 blocks.
!
         T11 = T( J1, J1 )
         T22 = T( J2, J2 )
!
!        Determine the transformation to perform the interchange.
!
         CALL SLARTG( T( J1, J2 ), T22-T11, CS, SN, TEMP )
!
!        Apply transformation to the matrix T.
!
         IF( J3<=N ) &
            CALL SROT( N-J1-1, T( J1, J3 ), LDT, T( J2, J3 ), LDT, CS, &
                       SN )
         CALL SROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
!
         T( J1, J1 ) = T22
         T( J2, J2 ) = T11
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL SROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
!
      ELSE
!
!        Swapping involves at least one 2-by-2 block.
!
!        Copy the diagonal block of order N1+N2 to the local array D
!        and compute its norm.
!
         ND = N1 + N2
         CALL SLACPY( 'Full', ND, ND, T( J1, J1 ), LDT, D, LDD )
         DNORM = SLANGE( 'Max', ND, ND, D, LDD, WORK )
!
!        Compute machine-dependent threshold for test for accepting
!        swap.
!
         EPS = SLAMCH( 'P' )
         SMLNUM = SLAMCH( 'S' ) / EPS
         THRESH = MAX( TEN*EPS*DNORM, SMLNUM )
!
!        Solve T11*X - X*T22 = scale*T12 for X.
!
         CALL SLASY2( .FALSE., .FALSE., -1, N1, N2, D, LDD, &
                      D( N1+1, N1+1 ), LDD, D( 1, N1+1 ), LDD, SCALE, X, &
                      LDX, XNORM, IERR )
!
!        Swap the adjacent diagonal blocks.
!
         K = N1 + N1 + N2 - 3
         GO TO ( 10, 20, 30 )K
!
   10    CONTINUE
!
!        N1 = 1, N2 = 2: generate elementary reflector H so that:
!
!        ( scale, X11, X12 ) H = ( 0, 0, * )
!
         U( 1 ) = SCALE
         U( 2 ) = X( 1, 1 )
         U( 3 ) = X( 1, 2 )
         CALL SLARFG( 3, U( 3 ), U, 1, TAU )
         U( 3 ) = ONE
         T11 = T( J1, J1 )
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL SLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL SLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
!
!        Test whether to reject swap.
!
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 3, &
             3 )-T11 ) )>THRESH )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL SLARFX( 'L', 3, N-J1+1, U, TAU, T( J1, J1 ), LDT, WORK )
         CALL SLARFX( 'R', J2, 3, U, TAU, T( 1, J1 ), LDT, WORK )
!
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J3, J3 ) = T11
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL SLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
!
   20    CONTINUE
!
!        N1 = 2, N2 = 1: generate elementary reflector H so that:
!
!        H (  -X11 ) = ( * )
!          (  -X21 ) = ( 0 )
!          ( scale ) = ( 0 )
!
         U( 1 ) = -X( 1, 1 )
         U( 2 ) = -X( 2, 1 )
         U( 3 ) = SCALE
         CALL SLARFG( 3, U( 1 ), U( 2 ), 1, TAU )
         U( 1 ) = ONE
         T33 = T( J3, J3 )
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL SLARFX( 'L', 3, 3, U, TAU, D, LDD, WORK )
         CALL SLARFX( 'R', 3, 3, U, TAU, D, LDD, WORK )
!
!        Test whether to reject swap.
!
         IF( MAX( ABS( D( 2, 1 ) ), ABS( D( 3, 1 ) ), ABS( D( 1, &
             1 )-T33 ) )>THRESH )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL SLARFX( 'R', J3, 3, U, TAU, T( 1, J1 ), LDT, WORK )
         CALL SLARFX( 'L', 3, N-J1, U, TAU, T( J1, J2 ), LDT, WORK )
!
         T( J1, J1 ) = T33
         T( J2, J1 ) = ZERO
         T( J3, J1 ) = ZERO
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL SLARFX( 'R', N, 3, U, TAU, Q( 1, J1 ), LDQ, WORK )
         END IF
         GO TO 40
!
   30    CONTINUE
!
!        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so
!        that:
!
!        H(2) H(1) (  -X11  -X12 ) = (  *  * )
!                  (  -X21  -X22 )   (  0  * )
!                  ( scale    0  )   (  0  0 )
!                  (    0  scale )   (  0  0 )
!
         U1( 1 ) = -X( 1, 1 )
         U1( 2 ) = -X( 2, 1 )
         U1( 3 ) = SCALE
         CALL SLARFG( 3, U1( 1 ), U1( 2 ), 1, TAU1 )
         U1( 1 ) = ONE
!
         TEMP = -TAU1*( X( 1, 2 )+U1( 2 )*X( 2, 2 ) )
         U2( 1 ) = -TEMP*U1( 2 ) - X( 2, 2 )
         U2( 2 ) = -TEMP*U1( 3 )
         U2( 3 ) = SCALE
         CALL SLARFG( 3, U2( 1 ), U2( 2 ), 1, TAU2 )
         U2( 1 ) = ONE
!
!        Perform swap provisionally on diagonal block in D.
!
         CALL SLARFX( 'L', 3, 4, U1, TAU1, D, LDD, WORK )
         CALL SLARFX( 'R', 4, 3, U1, TAU1, D, LDD, WORK )
         CALL SLARFX( 'L', 3, 4, U2, TAU2, D( 2, 1 ), LDD, WORK )
         CALL SLARFX( 'R', 4, 3, U2, TAU2, D( 1, 2 ), LDD, WORK )
!
!        Test whether to reject swap.
!
         IF( MAX( ABS( D( 3, 1 ) ), ABS( D( 3, 2 ) ), ABS( D( 4, 1 ) ), &
             ABS( D( 4, 2 ) ) )>THRESH )GO TO 50
!
!        Accept swap: apply transformation to the entire matrix T.
!
         CALL SLARFX( 'L', 3, N-J1+1, U1, TAU1, T( J1, J1 ), LDT, WORK )
         CALL SLARFX( 'R', J4, 3, U1, TAU1, T( 1, J1 ), LDT, WORK )
         CALL SLARFX( 'L', 3, N-J1+1, U2, TAU2, T( J2, J1 ), LDT, WORK )
         CALL SLARFX( 'R', J4, 3, U2, TAU2, T( 1, J2 ), LDT, WORK )
!
         T( J3, J1 ) = ZERO
         T( J3, J2 ) = ZERO
         T( J4, J1 ) = ZERO
         T( J4, J2 ) = ZERO
!
         IF( WANTQ ) THEN
!
!           Accumulate transformation in the matrix Q.
!
            CALL SLARFX( 'R', N, 3, U1, TAU1, Q( 1, J1 ), LDQ, WORK )
            CALL SLARFX( 'R', N, 3, U2, TAU2, Q( 1, J2 ), LDQ, WORK )
         END IF
!
   40    CONTINUE
!
         IF( N2==2 ) THEN
!
!           Standardize new 2-by-2 block T11
!
            CALL SLANV2( T( J1, J1 ), T( J1, J2 ), T( J2, J1 ), &
                         T( J2, J2 ), WR1, WI1, WR2, WI2, CS, SN )
            CALL SROT( N-J1-1, T( J1, J1+2 ), LDT, T( J2, J1+2 ), LDT, &
                       CS, SN )
            CALL SROT( J1-1, T( 1, J1 ), 1, T( 1, J2 ), 1, CS, SN )
            IF( WANTQ ) &
               CALL SROT( N, Q( 1, J1 ), 1, Q( 1, J2 ), 1, CS, SN )
         END IF
!
         IF( N1==2 ) THEN
!
!           Standardize new 2-by-2 block T22
!
            J3 = J1 + N2
            J4 = J3 + 1
            CALL SLANV2( T( J3, J3 ), T( J3, J4 ), T( J4, J3 ), &
                         T( J4, J4 ), WR1, WI1, WR2, WI2, CS, SN )
            IF( J3+2<=N ) &
               CALL SROT( N-J3-1, T( J3, J3+2 ), LDT, T( J4, J3+2 ), &
                          LDT, CS, SN )
            CALL SROT( J3-1, T( 1, J3 ), 1, T( 1, J4 ), 1, CS, SN )
            IF( WANTQ ) &
               CALL SROT( N, Q( 1, J3 ), 1, Q( 1, J4 ), 1, CS, SN )
         END IF
!
      END IF
      RETURN
!
!     Exit with INFO = 1 if swap was rejected.
!
   50 INFO = 1
      RETURN
!
!     End of SLAEXC
!
      END
      SUBROUTINE SLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, &
                         B, LDB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            LDB, LDX, N, NRHS
      REAL               ALPHA, BETA
!     ..
!     .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), &
                         X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  SLAGTM performs a matrix-vector product of the form
!
!     B := alpha * A * X + beta * B
!
!  where A is a tridiagonal matrix of order N, B and X are N by NRHS
!  matrices, and alpha and beta are real scalars, each of which may be
!  0., 1., or -1.
!
!  Arguments
!  =========
!
!  TRANS   (input) CHARACTER
!          Specifies the operation applied to A.
!          = 'N':  No transpose, B := alpha * A * X + beta * B
!          = 'T':  Transpose,    B := alpha * A'* X + beta * B
!          = 'C':  Conjugate transpose = Transpose
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices X and B.
!
!  ALPHA   (input) REAL
!          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,
!          it is assumed to be 0.
!
!  DL      (input) REAL array, dimension (N-1)
!          The (n-1) sub-diagonal elements of T.
!
!  D       (input) REAL array, dimension (N)
!          The diagonal elements of T.
!
!  DU      (input) REAL array, dimension (N-1)
!          The (n-1) super-diagonal elements of T.
!
!  X       (input) REAL array, dimension (LDX,NRHS)
!          The N by NRHS matrix X.
!  LDX     (input) INTEGER
!          The leading dimension of the array X.  LDX >= max(N,1).
!
!  BETA    (input) REAL
!          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,
!          it is assumed to be 1.
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the N by NRHS matrix B.
!          On exit, B is overwritten by the matrix expression
!          B := alpha * A * X + beta * B.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(N,1).
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
      IF( N==0 ) &
         RETURN
!
!     Multiply B by BETA if BETA/=1.
!
      IF( BETA==ZERO ) THEN
         DO 20 J = 1, NRHS
            DO 10 I = 1, N
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE IF( BETA==-ONE ) THEN
         DO 40 J = 1, NRHS
            DO 30 I = 1, N
               B( I, J ) = -B( I, J )
   30       CONTINUE
   40    CONTINUE
      END IF
!
      IF( ALPHA==ONE ) THEN
         IF( LSAME( TRANS, 'N' ) ) THEN
!
!           Compute B := B + A*X
!
            DO 60 J = 1, NRHS
               IF( N==1 ) THEN
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + &
                              DU( 1 )*X( 2, J )
                  B( N, J ) = B( N, J ) + DL( N-1 )*X( N-1, J ) + &
                              D( N )*X( N, J )
                  DO 50 I = 2, N - 1
                     B( I, J ) = B( I, J ) + DL( I-1 )*X( I-1, J ) + &
                                 D( I )*X( I, J ) + DU( I )*X( I+1, J )
   50             CONTINUE
               END IF
   60       CONTINUE
         ELSE
!
!           Compute B := B + A'*X
!
            DO 80 J = 1, NRHS
               IF( N==1 ) THEN
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + &
                              DL( 1 )*X( 2, J )
                  B( N, J ) = B( N, J ) + DU( N-1 )*X( N-1, J ) + &
                              D( N )*X( N, J )
                  DO 70 I = 2, N - 1
                     B( I, J ) = B( I, J ) + DU( I-1 )*X( I-1, J ) + &
                                 D( I )*X( I, J ) + DL( I )*X( I+1, J )
   70             CONTINUE
               END IF
   80       CONTINUE
         END IF
      ELSE IF( ALPHA==-ONE ) THEN
         IF( LSAME( TRANS, 'N' ) ) THEN
!
!           Compute B := B - A*X
!
            DO 100 J = 1, NRHS
               IF( N==1 ) THEN
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - &
                              DU( 1 )*X( 2, J )
                  B( N, J ) = B( N, J ) - DL( N-1 )*X( N-1, J ) - &
                              D( N )*X( N, J )
                  DO 90 I = 2, N - 1
                     B( I, J ) = B( I, J ) - DL( I-1 )*X( I-1, J ) - &
                                 D( I )*X( I, J ) - DU( I )*X( I+1, J )
   90             CONTINUE
               END IF
  100       CONTINUE
         ELSE
!
!           Compute B := B - A'*X
!
            DO 120 J = 1, NRHS
               IF( N==1 ) THEN
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               ELSE
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - &
                              DL( 1 )*X( 2, J )
                  B( N, J ) = B( N, J ) - DU( N-1 )*X( N-1, J ) - &
                              D( N )*X( N, J )
                  DO 110 I = 2, N - 1
                     B( I, J ) = B( I, J ) - DU( I-1 )*X( I-1, J ) - &
                                 D( I )*X( I, J ) - DL( I )*X( I+1, J )
  110             CONTINUE
               END IF
  120       CONTINUE
         END IF
      END IF
      RETURN
!
!     End of SLAGTM
!
      END
      SUBROUTINE SLAHQR( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI, &
                         ILOZ, IHIZ, Z, LDZ, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            WANTT, WANTZ
      INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, N
!     ..
!     .. Array Arguments ..
      REAL               H( LDH, * ), WI( * ), WR( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  SLAHQR is an auxiliary routine called by SHSEQR to update the
!  eigenvalues and Schur decomposition already computed by SHSEQR, by
!  dealing with the Hessenberg submatrix in rows and columns ILO to IHI.
!
!  Arguments
!  =========
!
!  WANTT   (input) LOGICAL
!          = .TRUE. : the full Schur form T is required;
!          = .FALSE.: only eigenvalues are required.
!
!  WANTZ   (input) LOGICAL
!          = .TRUE. : the matrix of Schur vectors Z is required;
!          = .FALSE.: Schur vectors are not required.
!
!  N       (input) INTEGER
!          The order of the matrix H.  N >= 0.
!
!  ILO     (input) INTEGER
!  IHI     (input) INTEGER
!          It is assumed that H is already upper quasi-triangular in
!          rows and columns IHI+1:N, and that H(ILO,ILO-1) = 0 (unless
!          ILO = 1). SLAHQR works primarily with the Hessenberg
!          submatrix in rows and columns ILO to IHI, but applies
!          transformations to all of H if WANTT is .TRUE..
!          1 <= ILO <= max(1,IHI); IHI <= N.
!
!  H       (input/output) REAL array, dimension (LDH,N)
!          On entry, the upper Hessenberg matrix H.
!          On exit, if WANTT is .TRUE., H is upper quasi-triangular in
!          rows and columns ILO:IHI, with any 2-by-2 diagonal blocks in
!          standard form. If WANTT is .FALSE., the contents of H are
!          unspecified on exit.
!
!  LDH     (input) INTEGER
!          The leading dimension of the array H. LDH >= max(1,N).
!
!  WR      (output) REAL array, dimension (N)
!  WI      (output) REAL array, dimension (N)
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
!  ILOZ    (input) INTEGER
!  IHIZ    (input) INTEGER
!          Specify the rows of Z to which transformations must be
!          applied if WANTZ is .TRUE..
!          1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!
!  Z       (input/output) REAL array, dimension (LDZ,N)
!          If WANTZ is .TRUE., on entry Z must contain the current
!          matrix Z of transformations accumulated by SHSEQR, and on
!          exit Z has been updated; transformations are applied only to
!          the submatrix Z(ILOZ:IHIZ,ILO:IHI).
!          If WANTZ is .FALSE., Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z. LDZ >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          > 0: SLAHQR failed to compute all the eigenvalues ILO to IHI
!               in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
!               elements i+1:ihi of WR and WI contain those eigenvalues
!               which have been successfully computed.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               DAT1, DAT2
      PARAMETER          ( DAT1 = 0.75E+0, DAT2 = -0.4375E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I1, I2, ITN, ITS, J, K, L, M, NH, NR, NZ
      REAL               CS, H00, H10, H11, H12, H21, H22, H33, H33S, &
                         H43H34, H44, H44S, OVFL, S, SMLNUM, SN, SUM, &
                         T1, T2, T3, TST1, ULP, UNFL, V1, V2, V3
!     ..
!     .. Local Arrays ..
      REAL               V( 3 ), WORK( 1 )
!     ..
!     .. External Functions ..
      REAL               SLAMCH, SLANHS
      EXTERNAL           SLAMCH, SLANHS
!     ..
!     .. External Subroutines ..
      EXTERNAL           SCOPY, SLABAD, SLANV2, SLARFG, SROT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Executable Statements ..
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N==0 ) &
         RETURN
      IF( ILO==IHI ) THEN
         WR( ILO ) = H( ILO, ILO )
         WI( ILO ) = ZERO
         RETURN
      END IF
!
      NH = IHI - ILO + 1
      NZ = IHIZ - ILOZ + 1
!
!     Set machine-dependent constants for the stopping criterion.
!     If norm(H) <= sqrt(OVFL), overflow should not occur.
!
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( NH / ULP )
!
!     I1 and I2 are the indices of the first row and last column of H
!     to which transformations must be applied. If eigenvalues only are
!     being computed, I1 and I2 are set inside the main loop.
!
      IF( WANTT ) THEN
         I1 = 1
         I2 = N
      END IF
!
!     ITN is the total number of QR iterations allowed.
!
      ITN = 30*NH
!
!     The main loop begins here. I is the loop index and decreases from
!     IHI to ILO in steps of 1 or 2. Each iteration of the loop works
!     with the active submatrix in rows and columns L to I.
!     Eigenvalues I+1 to IHI have already converged. Either L = ILO or
!     H(L,L-1) is negligible so that the matrix splits.
!
      I = IHI
   10 CONTINUE
      L = ILO
      IF( I<ILO ) &
         GO TO 150
!
!     Perform QR iterations on rows and columns ILO to I until a
!     submatrix of order 1 or 2 splits off at the bottom because a
!     subdiagonal element has become negligible.
!
      DO 130 ITS = 0, ITN
!
!        Look for a single small subdiagonal element.
!
         DO 20 K = I, L + 1, -1
            TST1 = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
            IF( TST1==ZERO ) &
               TST1 = SLANHS( '1', I-L+1, H( L, L ), LDH, WORK )
            IF( ABS( H( K, K-1 ) )<=MAX( ULP*TST1, SMLNUM ) ) &
               GO TO 30
   20    CONTINUE
   30    CONTINUE
         L = K
         IF( L>ILO ) THEN
!
!           H(L,L-1) is negligible
!
            H( L, L-1 ) = ZERO
         END IF
!
!        Exit from loop if a submatrix of order 1 or 2 has split off.
!
         IF( L>=I-1 ) &
            GO TO 140
!
!        Now the active submatrix is in rows and columns L to I. If
!        eigenvalues only are being computed, only the active submatrix
!        need be transformed.
!
         IF( .NOT.WANTT ) THEN
            I1 = L
            I2 = I
         END IF
!
         IF( ITS==10 .OR. ITS==20 ) THEN
!
!           Exceptional shift.
!
            S = ABS( H( I, I-1 ) ) + ABS( H( I-1, I-2 ) )
            H44 = DAT1*S
            H33 = H44
            H43H34 = DAT2*S*S
         ELSE
!
!           Prepare to use Wilkinson's double shift
!
            H44 = H( I, I )
            H33 = H( I-1, I-1 )
            H43H34 = H( I, I-1 )*H( I-1, I )
         END IF
!
!        Look for two consecutive small subdiagonal elements.
!
         DO 40 M = I - 2, L, -1
!
!           Determine the effect of starting the double-shift QR
!           iteration at row M, and see if this would make H(M,M-1)
!           negligible.
!
            H11 = H( M, M )
            H22 = H( M+1, M+1 )
            H21 = H( M+1, M )
            H12 = H( M, M+1 )
            H44S = H44 - H11
            H33S = H33 - H11
            V1 = ( H33S*H44S-H43H34 ) / H21 + H12
            V2 = H22 - H11 - H33S - H44S
            V3 = H( M+2, M+1 )
            S = ABS( V1 ) + ABS( V2 ) + ABS( V3 )
            V1 = V1 / S
            V2 = V2 / S
            V3 = V3 / S
            V( 1 ) = V1
            V( 2 ) = V2
            V( 3 ) = V3
            IF( M==L ) &
               GO TO 50
            H00 = H( M-1, M-1 )
            H10 = H( M, M-1 )
            TST1 = ABS( V1 )*( ABS( H00 )+ABS( H11 )+ABS( H22 ) )
            IF( ABS( H10 )*( ABS( V2 )+ABS( V3 ) )<=ULP*TST1 ) &
               GO TO 50
   40    CONTINUE
   50    CONTINUE
!
!        Double-shift QR step
!
         DO 120 K = M, I - 1
!
!           The first iteration of this loop determines a reflection G
!           from the vector V and applies it from left and right to H,
!           thus creating a nonzero bulge below the subdiagonal.
!
!           Each subsequent iteration determines a reflection G to
!           restore the Hessenberg form in the (K-1)th column, and thus
!           chases the bulge one step toward the bottom of the active
!           submatrix. NR is the order of G.
!
            NR = MIN( 3, I-K+1 )
            IF( K>M ) &
               CALL SCOPY( NR, H( K, K-1 ), 1, V, 1 )
            CALL SLARFG( NR, V( 1 ), V( 2 ), 1, T1 )
            IF( K>M ) THEN
               H( K, K-1 ) = V( 1 )
               H( K+1, K-1 ) = ZERO
               IF( K<I-1 ) &
                  H( K+2, K-1 ) = ZERO
            ELSE IF( M>L ) THEN
               H( K, K-1 ) = -H( K, K-1 )
            END IF
            V2 = V( 2 )
            T2 = T1*V2
            IF( NR==3 ) THEN
               V3 = V( 3 )
               T3 = T1*V3
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 60 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J ) + V3*H( K+2, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
                  H( K+2, J ) = H( K+2, J ) - SUM*T3
   60          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 70 J = I1, MIN( K+3, I )
                  SUM = H( J, K ) + V2*H( J, K+1 ) + V3*H( J, K+2 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
                  H( J, K+2 ) = H( J, K+2 ) - SUM*T3
   70          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 80 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 ) + V3*Z( J, K+2 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
                     Z( J, K+2 ) = Z( J, K+2 ) - SUM*T3
   80             CONTINUE
               END IF
            ELSE IF( NR==2 ) THEN
!
!              Apply G from the left to transform the rows of the matrix
!              in columns K to I2.
!
               DO 90 J = K, I2
                  SUM = H( K, J ) + V2*H( K+1, J )
                  H( K, J ) = H( K, J ) - SUM*T1
                  H( K+1, J ) = H( K+1, J ) - SUM*T2
   90          CONTINUE
!
!              Apply G from the right to transform the columns of the
!              matrix in rows I1 to min(K+3,I).
!
               DO 100 J = I1, I
                  SUM = H( J, K ) + V2*H( J, K+1 )
                  H( J, K ) = H( J, K ) - SUM*T1
                  H( J, K+1 ) = H( J, K+1 ) - SUM*T2
  100          CONTINUE
!
               IF( WANTZ ) THEN
!
!                 Accumulate transformations in the matrix Z
!
                  DO 110 J = ILOZ, IHIZ
                     SUM = Z( J, K ) + V2*Z( J, K+1 )
                     Z( J, K ) = Z( J, K ) - SUM*T1
                     Z( J, K+1 ) = Z( J, K+1 ) - SUM*T2
  110             CONTINUE
               END IF
            END IF
  120    CONTINUE
!
  130 CONTINUE
!
!     Failure to converge in remaining number of iterations
!
      INFO = I
      RETURN
!
  140 CONTINUE
!
      IF( L==I ) THEN
!
!        H(I,I-1) is negligible: one eigenvalue has converged.
!
         WR( I ) = H( I, I )
         WI( I ) = ZERO
      ELSE IF( L==I-1 ) THEN
!
!        H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
!
!        Transform the 2-by-2 submatrix to standard Schur form,
!        and compute and store the eigenvalues.
!
         CALL SLANV2( H( I-1, I-1 ), H( I-1, I ), H( I, I-1 ), &
                      H( I, I ), WR( I-1 ), WI( I-1 ), WR( I ), WI( I ), &
                      CS, SN )
!
         IF( WANTT ) THEN
!
!           Apply the transformation to the rest of H.
!
            IF( I2>I ) &
               CALL SROT( I2-I, H( I-1, I+1 ), LDH, H( I, I+1 ), LDH, &
                          CS, SN )
            CALL SROT( I-I1-1, H( I1, I-1 ), 1, H( I1, I ), 1, CS, SN )
         END IF
         IF( WANTZ ) THEN
!
!           Apply the transformation to Z.
!
            CALL SROT( NZ, Z( ILOZ, I-1 ), 1, Z( ILOZ, I ), 1, CS, SN )
         END IF
      END IF
!
!     Decrement number of remaining iterations, and return to start of
!     the main loop with new value of I.
!
      ITN = ITN - ITS
      I = L - 1
      GO TO 10
!
  150 CONTINUE
      RETURN
!
!     End of SLAHQR
!
      END
      SUBROUTINE SLALN2( LTRANS, NA, NW, SMIN, CA, A, LDA, D1, D2, B, &
                         LDB, WR, WI, X, LDX, SCALE, XNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            LTRANS
      INTEGER            INFO, LDA, LDB, LDX, NA, NW
      REAL               CA, D1, D2, SCALE, SMIN, WI, WR, XNORM
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  SLALN2 solves a system of the form  (ca A - w D ) X = s B
!  or (ca A' - w D) X = s B   with possible scaling ("s") and
!  perturbation of A.  (A' means A-transpose.)
!
!  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
!  real diagonal matrix, w is a real or complex value, and X and B are
!  NA x 1 matrices -- real if w is real, complex if w is complex.  NA
!  may be 1 or 2.
!
!  If w is complex, X and B are represented as NA x 2 matrices,
!  the first column of each being the real part and the second
!  being the imaginary part.
!
!  "s" is a scaling factor (<= 1), computed by SLALN2, which is
!  so chosen that X can be computed without overflow.  X is further
!  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
!  than overflow.
!
!  If both singular values of (ca A - w D) are less than SMIN,
!  SMIN*identity will be used instead of (ca A - w D).  If only one
!  singular value is less than SMIN, one element of (ca A - w D) will be
!  perturbed enough to make the smallest singular value roughly SMIN.
!  If both singular values are at least SMIN, (ca A - w D) will not be
!  perturbed.  In any case, the perturbation will be at most some small
!  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
!  are computed by infinity-norm approximations, and thus will only be
!  correct to a factor of 2 or so.
!
!  Note: all input quantities are assumed to be smaller than overflow
!  by a reasonable factor.  (See BIGNUM.)
!
!  Arguments
!  ==========
!
!  LTRANS  (input) LOGICAL
!          =.TRUE.:  A-transpose will be used.
!          =.FALSE.: A will be used (not transposed.)
!
!  NA      (input) INTEGER
!          The size of the matrix A.  It may (only) be 1 or 2.
!
!  NW      (input) INTEGER
!          1 if "w" is real, 2 if "w" is complex.  It may only be 1
!          or 2.
!
!  SMIN    (input) REAL
!          The desired lower bound on the singular values of A.  This
!          should be a safe distance away from underflow or overflow,
!          say, between (underflow/machine precision) and  (machine
!          precision * overflow ).  (See BIGNUM and ULP.)
!
!  CA      (input) REAL
!          The coefficient c, which A is multiplied by.
!
!  A       (input) REAL array, dimension (LDA,NA)
!          The NA x NA matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of A.  It must be at least NA.
!
!  D1      (input) REAL
!          The 1,1 element in the diagonal matrix D.
!
!  D2      (input) REAL
!          The 2,2 element in the diagonal matrix D.  Not used if NW=1.
!
!  B       (input) REAL array, dimension (LDB,NW)
!          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
!          complex), column 1 contains the real part of B and column 2
!          contains the imaginary part.
!
!  LDB     (input) INTEGER
!          The leading dimension of B.  It must be at least NA.
!
!  WR      (input) REAL
!          The real part of the scalar "w".
!
!  WI      (input) REAL
!          The imaginary part of the scalar "w".  Not used if NW=1.
!
!  X       (output) REAL array, dimension (LDX,NW)
!          The NA x NW matrix X (unknowns), as computed by SLALN2.
!          If NW=2 ("w" is complex), on exit, column 1 will contain
!          the real part of X and column 2 will contain the imaginary
!          part.
!
!  LDX     (input) INTEGER
!          The leading dimension of X.  It must be at least NA.
!
!  SCALE   (output) REAL
!          The scale factor that B must be multiplied by to insure
!          that overflow does not occur when computing X.  Thus,
!          (ca A - w D) X  will be SCALE*B, not B (ignoring
!          perturbations of A.)  It will be at most 1.
!
!  XNORM   (output) REAL
!          The infinity-norm of X, when X is regarded as an NA x NW
!          real matrix.
!
!  INFO    (output) INTEGER
!          An error flag.  It will be set to zero if no error occurs,
!          a negative number if an argument is in error, or a positive
!          number if  ca A - w D  had to be perturbed.
!          The possible values are:
!          = 0: No error occurred, and (ca A - w D) did not have to be
!                 perturbed.
!          = 1: (ca A - w D) had to be perturbed to make its smallest
!               (or only) singular value greater than SMIN.
!          NOTE: In the interests of speed, this routine does not
!                check the inputs for errors.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
!     ..
!     .. Local Scalars ..
      INTEGER            ICMAX, J
      REAL               BBND, BI1, BI2, BIGNUM, BNORM, BR1, BR2, CI21, &
                         CI22, CMAX, CNORM, CR21, CR22, CSI, CSR, LI21, &
                         LR21, SMINI, SMLNUM, TEMP, U22ABS, UI11, UI11R, &
                         UI12, UI12S, UI22, UR11, UR11R, UR12, UR12S, &
                         UR22, XI1, XI2, XR1, XR2
!     ..
!     .. Local Arrays ..
      LOGICAL            CSWAP( 4 ), RSWAP( 4 )
      INTEGER            IPIVOT( 4, 4 )
      REAL               CI( 2, 2 ), CIV( 4 ), CR( 2, 2 ), CRV( 4 )
!     ..
!     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Equivalences ..
      EQUIVALENCE        ( CI( 1, 1 ), CIV( 1 ) ), &
                         ( CR( 1, 1 ), CRV( 1 ) )
!     ..
!     .. Data statements ..
      DATA               CSWAP / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               RSWAP / .FALSE., .TRUE., .FALSE., .TRUE. /
      DATA               IPIVOT / 1, 2, 3, 4, 2, 1, 4, 3, 3, 4, 1, 2, 4, &
                         3, 2, 1 /
!     ..
!     .. Executable Statements ..
!
!     Compute BIGNUM
!
      SMLNUM = TWO*SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      SMINI = MAX( SMIN, SMLNUM )
!
!     Don't check for input errors
!
      INFO = 0
!
!     Standard Initializations
!
      SCALE = ONE
!
      IF( NA==1 ) THEN
!
!        1 x 1  (i.e., scalar) system   C X = B
!
         IF( NW==1 ) THEN
!
!           Real 1x1 system.
!
!           C = ca A - w D
!
            CSR = CA*A( 1, 1 ) - WR*D1
            CNORM = ABS( CSR )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( CNORM<SMINI ) THEN
               CSR = SMINI
               CNORM = SMINI
               INFO = 1
            END IF
!
!           Check scaling for  X = B / C
!
            BNORM = ABS( B( 1, 1 ) )
            IF( CNORM<ONE .AND. BNORM>ONE ) THEN
               IF( BNORM>BIGNUM*CNORM ) &
                  SCALE = ONE / BNORM
            END IF
!
!           Compute X
!
            X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / CSR
            XNORM = ABS( X( 1, 1 ) )
         ELSE
!
!           Complex 1x1 system (w is complex)
!
!           C = ca A - w D
!
            CSR = CA*A( 1, 1 ) - WR*D1
            CSI = -WI*D1
            CNORM = ABS( CSR ) + ABS( CSI )
!
!           If | C | < SMINI, use C = SMINI
!
            IF( CNORM<SMINI ) THEN
               CSR = SMINI
               CSI = ZERO
               CNORM = SMINI
               INFO = 1
            END IF
!
!           Check scaling for  X = B / C
!
            BNORM = ABS( B( 1, 1 ) ) + ABS( B( 1, 2 ) )
            IF( CNORM<ONE .AND. BNORM>ONE ) THEN
               IF( BNORM>BIGNUM*CNORM ) &
                  SCALE = ONE / BNORM
            END IF
!
!           Compute X
!
            CALL SLADIV( SCALE*B( 1, 1 ), SCALE*B( 1, 2 ), CSR, CSI, &
                         X( 1, 1 ), X( 1, 2 ) )
            XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
         END IF
!
      ELSE
!
!        2x2 System
!
!        Compute the real part of  C = ca A - w D  (or  ca A' - w D )
!
         CR( 1, 1 ) = CA*A( 1, 1 ) - WR*D1
         CR( 2, 2 ) = CA*A( 2, 2 ) - WR*D2
         IF( LTRANS ) THEN
            CR( 1, 2 ) = CA*A( 2, 1 )
            CR( 2, 1 ) = CA*A( 1, 2 )
         ELSE
            CR( 2, 1 ) = CA*A( 2, 1 )
            CR( 1, 2 ) = CA*A( 1, 2 )
         END IF
!
         IF( NW==1 ) THEN
!
!           Real 2x2 system  (w is real)
!
!           Find the largest element in C
!
            CMAX = ZERO
            ICMAX = 0
!
            DO 10 J = 1, 4
               IF( ABS( CRV( J ) )>CMAX ) THEN
                  CMAX = ABS( CRV( J ) )
                  ICMAX = J
               END IF
   10       CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( CMAX<SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) ), ABS( B( 2, 1 ) ) )
               IF( SMINI<ONE .AND. BNORM>ONE ) THEN
                  IF( BNORM>BIGNUM*SMINI ) &
                     SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            UR11 = CRV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            UR11R = ONE / UR11
            LR21 = UR11R*CR21
            UR22 = CR22 - UR12*LR21
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( ABS( UR22 )<SMINI ) THEN
               UR22 = SMINI
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR1 = B( 2, 1 )
               BR2 = B( 1, 1 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
            END IF
            BR2 = BR2 - LR21*BR1
            BBND = MAX( ABS( BR1*( UR22*UR11R ) ), ABS( BR2 ) )
            IF( BBND>ONE .AND. ABS( UR22 )<ONE ) THEN
               IF( BBND>=BIGNUM*ABS( UR22 ) ) &
                  SCALE = ONE / BBND
            END IF
!
            XR2 = ( BR2*SCALE ) / UR22
            XR1 = ( SCALE*BR1 )*UR11R - XR2*( UR11R*UR12 )
            IF( CSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
            END IF
            XNORM = MAX( ABS( XR1 ), ABS( XR2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( XNORM>ONE .AND. CMAX>ONE ) THEN
               IF( XNORM>BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         ELSE
!
!           Complex 2x2 system  (w is complex)
!
!           Find the largest element in C
!
            CI( 1, 1 ) = -WI*D1
            CI( 2, 1 ) = ZERO
            CI( 1, 2 ) = ZERO
            CI( 2, 2 ) = -WI*D2
            CMAX = ZERO
            ICMAX = 0
!
            DO 20 J = 1, 4
               IF( ABS( CRV( J ) )+ABS( CIV( J ) )>CMAX ) THEN
                  CMAX = ABS( CRV( J ) ) + ABS( CIV( J ) )
                  ICMAX = J
               END IF
   20       CONTINUE
!
!           If norm(C) < SMINI, use SMINI*identity.
!
            IF( CMAX<SMINI ) THEN
               BNORM = MAX( ABS( B( 1, 1 ) )+ABS( B( 1, 2 ) ), &
                       ABS( B( 2, 1 ) )+ABS( B( 2, 2 ) ) )
               IF( SMINI<ONE .AND. BNORM>ONE ) THEN
                  IF( BNORM>BIGNUM*SMINI ) &
                     SCALE = ONE / BNORM
               END IF
               TEMP = SCALE / SMINI
               X( 1, 1 ) = TEMP*B( 1, 1 )
               X( 2, 1 ) = TEMP*B( 2, 1 )
               X( 1, 2 ) = TEMP*B( 1, 2 )
               X( 2, 2 ) = TEMP*B( 2, 2 )
               XNORM = TEMP*BNORM
               INFO = 1
               RETURN
            END IF
!
!           Gaussian elimination with complete pivoting.
!
            UR11 = CRV( ICMAX )
            UI11 = CIV( ICMAX )
            CR21 = CRV( IPIVOT( 2, ICMAX ) )
            CI21 = CIV( IPIVOT( 2, ICMAX ) )
            UR12 = CRV( IPIVOT( 3, ICMAX ) )
            UI12 = CIV( IPIVOT( 3, ICMAX ) )
            CR22 = CRV( IPIVOT( 4, ICMAX ) )
            CI22 = CIV( IPIVOT( 4, ICMAX ) )
            IF( ICMAX==1 .OR. ICMAX==4 ) THEN
!
!              Code when off-diagonals of pivoted C are real
!
               IF( ABS( UR11 )>ABS( UI11 ) ) THEN
                  TEMP = UI11 / UR11
                  UR11R = ONE / ( UR11*( ONE+TEMP**2 ) )
                  UI11R = -TEMP*UR11R
               ELSE
                  TEMP = UR11 / UI11
                  UI11R = -ONE / ( UI11*( ONE+TEMP**2 ) )
                  UR11R = -TEMP*UI11R
               END IF
               LR21 = CR21*UR11R
               LI21 = CR21*UI11R
               UR12S = UR12*UR11R
               UI12S = UR12*UI11R
               UR22 = CR22 - UR12*LR21
               UI22 = CI22 - UR12*LI21
            ELSE
!
!              Code when diagonals of pivoted C are real
!
               UR11R = ONE / UR11
               UI11R = ZERO
               LR21 = CR21*UR11R
               LI21 = CI21*UR11R
               UR12S = UR12*UR11R
               UI12S = UI12*UR11R
               UR22 = CR22 - UR12*LR21 + UI12*LI21
               UI22 = -UR12*LI21 - UI12*LR21
            END IF
            U22ABS = ABS( UR22 ) + ABS( UI22 )
!
!           If smaller pivot < SMINI, use SMINI
!
            IF( U22ABS<SMINI ) THEN
               UR22 = SMINI
               UI22 = ZERO
               INFO = 1
            END IF
            IF( RSWAP( ICMAX ) ) THEN
               BR2 = B( 1, 1 )
               BR1 = B( 2, 1 )
               BI2 = B( 1, 2 )
               BI1 = B( 2, 2 )
            ELSE
               BR1 = B( 1, 1 )
               BR2 = B( 2, 1 )
               BI1 = B( 1, 2 )
               BI2 = B( 2, 2 )
            END IF
            BR2 = BR2 - LR21*BR1 + LI21*BI1
            BI2 = BI2 - LI21*BR1 - LR21*BI1
            BBND = MAX( ( ABS( BR1 )+ABS( BI1 ) )* &
                   ( U22ABS*( ABS( UR11R )+ABS( UI11R ) ) ), &
                   ABS( BR2 )+ABS( BI2 ) )
            IF( BBND>ONE .AND. U22ABS<ONE ) THEN
               IF( BBND>=BIGNUM*U22ABS ) THEN
                  SCALE = ONE / BBND
                  BR1 = SCALE*BR1
                  BI1 = SCALE*BI1
                  BR2 = SCALE*BR2
                  BI2 = SCALE*BI2
               END IF
            END IF
!
            CALL SLADIV( BR2, BI2, UR22, UI22, XR2, XI2 )
            XR1 = UR11R*BR1 - UI11R*BI1 - UR12S*XR2 + UI12S*XI2
            XI1 = UI11R*BR1 + UR11R*BI1 - UI12S*XR2 - UR12S*XI2
            IF( CSWAP( ICMAX ) ) THEN
               X( 1, 1 ) = XR2
               X( 2, 1 ) = XR1
               X( 1, 2 ) = XI2
               X( 2, 2 ) = XI1
            ELSE
               X( 1, 1 ) = XR1
               X( 2, 1 ) = XR2
               X( 1, 2 ) = XI1
               X( 2, 2 ) = XI2
            END IF
            XNORM = MAX( ABS( XR1 )+ABS( XI1 ), ABS( XR2 )+ABS( XI2 ) )
!
!           Further scaling if  norm(A) norm(X) > overflow
!
            IF( XNORM>ONE .AND. CMAX>ONE ) THEN
               IF( XNORM>BIGNUM / CMAX ) THEN
                  TEMP = CMAX / BIGNUM
                  X( 1, 1 ) = TEMP*X( 1, 1 )
                  X( 2, 1 ) = TEMP*X( 2, 1 )
                  X( 1, 2 ) = TEMP*X( 1, 2 )
                  X( 2, 2 ) = TEMP*X( 2, 2 )
                  XNORM = TEMP*XNORM
                  SCALE = TEMP*SCALE
               END IF
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of SLALN2
!
      END
      REAL             FUNCTION SLAMCH( CMACH )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          CMACH
!     ..
!
!  Purpose
!  =======
!
!  SLAMCH determines single precision machine parameters.
!
!  Arguments
!  =========
!
!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by SLAMCH:
!          = 'E' or 'e',   SLAMCH := eps
!          = 'S' or 's ,   SLAMCH := sfmin
!          = 'B' or 'b',   SLAMCH := base
!          = 'P' or 'p',   SLAMCH := eps*base
!          = 'N' or 'n',   SLAMCH := t
!          = 'R' or 'r',   SLAMCH := rnd
!          = 'M' or 'm',   SLAMCH := emin
!          = 'U' or 'u',   SLAMCH := rmin
!          = 'L' or 'l',   SLAMCH := emax
!          = 'O' or 'o',   SLAMCH := rmax
!
!          where
!
!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      REAL               BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN, &
                         RND, SFMIN, SMALL, T
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAMC2
!     ..
!     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN, &
                         EMAX, RMAX, PREC
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL SLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL>=SFMIN ) THEN
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
!
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
!
      SLAMCH = RMACH
      RETURN
      END
      SUBROUTINE SLAMC1( BETA, T, RND, IEEE1 )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
!     ..
!
!  Purpose
!  =======
!
!  SLAMC1 determines the machine parameters given by BETA, T, RND, and
!  IEEE1.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  IEEE1   (output) LOGICAL
!          Specifies whether rounding appears to be done in the IEEE
!          'round to nearest' style.
!
!  Further Details
!  ===============
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
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      REAL               A, B, C, F, ONE, QTR, SAVEC, T1, T2
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         A = 1
         C = 1
!
!+       WHILE( C==ONE )LOOP
   10    CONTINUE
         IF( C==ONE ) THEN
            A = 2*A
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 10
         END IF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
         B = 1
         C = SLAMC3( A, B )
!
!+       WHILE( C==A )LOOP
   20    CONTINUE
         IF( C==A ) THEN
            B = 2*B
            C = SLAMC3( A, B )
            GO TO 20
         END IF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
         QTR = ONE / 4
         SAVEC = C
         C = SLAMC3( C, -A )
         LBETA = C + QTR
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
         B = LBETA
         F = SLAMC3( B / 2, -B / 100 )
         C = SLAMC3( F, A )
         IF( C==A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = SLAMC3( B / 2, B / 100 )
         C = SLAMC3( F, A )
         IF( ( LRND ) .AND. ( C==A ) ) &
            LRND = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
         T1 = SLAMC3( B / 2, A )
         T2 = SLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1==A ) .AND. ( T2>SAVEC ) .AND. LRND
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
         LT = 0
         A = 1
         C = 1
!
!+       WHILE( C==ONE )LOOP
   30    CONTINUE
         IF( C==ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = SLAMC3( A, ONE )
            C = SLAMC3( C, -A )
            GO TO 30
         END IF
!+       END WHILE
!
      END IF
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
!
!     End of SLAMC1
!
      END
      SUBROUTINE SLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      REAL               EPS, RMAX, RMIN
!     ..
!
!  Purpose
!  =======
!
!  SLAMC2 determines the machine parameters specified in its argument
!  list.
!
!  Arguments
!  =========
!
!  BETA    (output) INTEGER
!          The base of the machine.
!
!  T       (output) INTEGER
!          The number of ( BETA ) digits in the mantissa.
!
!  RND     (output) LOGICAL
!          Specifies whether proper rounding  ( RND = .TRUE. )  or
!          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!          be a reliable guide to the way in which the machine performs
!          its arithmetic.
!
!  EPS     (output) REAL
!          The smallest positive number such that
!
!             fl( 1.0 - EPS ) < 1.0,
!
!          where fl denotes the computed value.
!
!  EMIN    (output) INTEGER
!          The minimum exponent before (gradual) underflow occurs.
!
!  RMIN    (output) REAL
!          The smallest normalized number for the machine, given by
!          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!          of BETA.
!
!  EMAX    (output) INTEGER
!          The maximum exponent before overflow occurs.
!
!  RMAX    (output) REAL
!          The largest positive number for the machine, given by
!          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!          value of BETA.
!
!  Further Details
!  ===============
!
!  The computation of  EPS  is based on a routine PARANOIA by
!  W. Kahan of the University of California at Berkeley.
!
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT, &
                         NGNMIN, NGPMIN
      REAL               A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE, &
                         SIXTH, SMALL, THIRD, TWO, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAMC1, SLAMC4, SLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX, &
                         LRMIN, LT
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  SLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL SLAMC1( LBETA, LT, LRND, LIEEE1 )
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
         SIXTH = SLAMC3( B, -HALF )
         THIRD = SLAMC3( SIXTH, SIXTH )
         B = SLAMC3( THIRD, -HALF )
         B = SLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B<LEPS ) &
            B = LEPS
!
         LEPS = 1
!
!+       WHILE( ( LEPS>B ).AND.( B>ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS>B ) .AND. ( B>ZERO ) ) THEN
            LEPS = B
            C = SLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = SLAMC3( HALF, -C )
            B = SLAMC3( HALF, C )
            C = SLAMC3( HALF, -B )
            B = SLAMC3( HALF, C )
            GO TO 10
         END IF
!+       END WHILE
!
         IF( A<LEPS ) &
            LEPS = A
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = SLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = SLAMC3( ONE, SMALL )
         CALL SLAMC4( NGPMIN, ONE, LBETA )
         CALL SLAMC4( NGNMIN, -ONE, LBETA )
         CALL SLAMC4( GPMIN, A, LBETA )
         CALL SLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
!
         IF( ( NGPMIN==NGNMIN ) .AND. ( GPMIN==GNMIN ) ) THEN
            IF( NGPMIN==GPMIN ) THEN
               LEMIN = NGPMIN
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN )==3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( NGPMIN==GPMIN ) .AND. ( NGNMIN==GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN )==1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE IF( ( ABS( NGPMIN-NGNMIN )==1 ) .AND. &
                  ( GPMIN==GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) )==3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
!            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
!
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
!         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
!**
! Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine SLAMC1. A true IEEE machine should have both  things
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
            LRMIN = SLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
!
!        Finally, call SLAMC5 to compute EMAX and RMAX.
!
         CALL SLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
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
!
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-', &
            '  EMIN = ', I8, / &
            ' If, after inspection, the value EMIN looks', &
            ' acceptable please comment out ', &
            / ' the IF block as marked within the code of routine', &
            ' SLAMC2,', / ' otherwise supply EMIN explicitly.', / )
!
!     End of SLAMC2
!
      END
      REAL             FUNCTION SLAMC3( A, B )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               A, B
!     ..
!
!  Purpose
!  =======
!
!  SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!  the addition of  A  and  B ,  for use in situations where optimizers
!  might hold one of these in a register.
!
!  Arguments
!  =========
!
!  A, B    (input) REAL
!          The values A and B.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      SLAMC3 = A + B
!
      RETURN
!
!     End of SLAMC3
!
      END
      SUBROUTINE SLAMC4( EMIN, START, BASE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      REAL               START
!     ..
!
!  Purpose
!  =======
!
!  SLAMC4 is a service routine for SLAMC2.
!
!  Arguments
!  =========
!
!  EMIN    (output) EMIN
!          The minimum exponent before (gradual) underflow, computed by
!          setting A = START and dividing by BASE until the previous A
!          can not be recovered.
!
!  START   (input) REAL
!          The starting point for determining EMIN.
!
!  BASE    (input) INTEGER
!          The base of the machine.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I
      REAL               A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Executable Statements ..
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = SLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!+    WHILE( ( C1==A ).AND.( C2==A ).AND.
!    $       ( D1==A ).AND.( D2==A )      )LOOP
   10 CONTINUE
      IF( ( C1==A ) .AND. ( C2==A ) .AND. ( D1==A ) .AND. &
          ( D2==A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = SLAMC3( A / BASE, ZERO )
         C1 = SLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = SLAMC3( A*RBASE, ZERO )
         C2 = SLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
!+    END WHILE
!
      RETURN
!
!     End of SLAMC4
!
      END
      SUBROUTINE SLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      REAL               RMAX
!     ..
!
!  Purpose
!  =======
!
!  SLAMC5 attempts to compute RMAX, the largest machine floating-point
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!  approximately to a power of 2.  It will fail on machines where this
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is
!  too large (i.e. too close to zero), probably with overflow.
!
!  Arguments
!  =========
!
!  BETA    (input) INTEGER
!          The base of floating-point arithmetic.
!
!  P       (input) INTEGER
!          The number of base BETA digits in the mantissa of a
!          floating-point value.
!
!  EMIN    (input) INTEGER
!          The minimum exponent before (gradual) underflow.
!
!  IEEE    (input) LOGICAL
!          A logical flag specifying whether or not the arithmetic
!          system is thought to comply with the IEEE standard.
!
!  EMAX    (output) INTEGER
!          The largest exponent before overflow
!
!  RMAX    (output) REAL
!          The largest machine floating-point number.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      REAL               OLDY, RECBAS, Y, Z
!     ..
!     .. External Functions ..
      REAL               SLAMC3
      EXTERNAL           SLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD
!     ..
!     .. Executable Statements ..
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
      IF( TRY<=( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP==-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
      IF( ( UEXP+EMIN )>( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
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
      IF( ( MOD( NBITS, 2 )==1 ) .AND. ( BETA==2 ) ) THEN
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
      IF( IEEE ) THEN
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         EMAX = EMAX - 1
      END IF
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y<ONE ) &
            OLDY = Y
         Y = SLAMC3( Y, Z )
   20 CONTINUE
      IF( Y>=ONE ) &
         Y = OLDY
!
!     Now multiply by BETA**EMAX to get RMAX.
!
      DO 30 I = 1, EMAX
         Y = SLAMC3( Y*BETA, ZERO )
   30 CONTINUE
!
      RMAX = Y
      RETURN
!
!     End of SLAMC5
!
      END
      REAL             FUNCTION SLANGE( NORM, M, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLANGE  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real matrix A.
!
!  Description
!  ===========
!
!  SLANGE returns the value
!
!     SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in SLANGE as described
!          above.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.  When M = 0,
!          SLANGE is set to zero.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.  When N = 0,
!          SLANGE is set to zero.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The m by n matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(M,1).
!
!  WORK    (workspace) REAL array, dimension (LWORK),
!          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      REAL               SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( MIN( M, N )==0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM=='1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL SLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      SLANGE = VALUE
      RETURN
!
!     End of SLANGE
!
      END
      REAL             FUNCTION SLANHS( NORM, N, A, LDA, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLANHS  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  Hessenberg matrix A.
!
!  Description
!  ===========
!
!  SLANHS returns the value
!
!     SLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in SLANHS as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, SLANHS is
!          set to zero.
!
!  A       (input) REAL array, dimension (LDA,N)
!          The n by n upper Hessenberg matrix A; the part of A below the
!          first sub-diagonal is not referenced.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).
!
!  WORK    (workspace) REAL array, dimension (LWORK),
!          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
!          referenced.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
      REAL               SCALE, SUM, VALUE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLASSQ
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N==0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, MIN( N, J+1 )
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM=='1' ) ) THEN
!
!        Find norm1(A).
!
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
!
!        Find normI(A).
!
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL SLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
!
      SLANHS = VALUE
      RETURN
!
!     End of SLANHS
!
      END
      REAL             FUNCTION SLANST( NORM, N, D, E )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
!     ..
!     .. Array Arguments ..
      REAL               D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  SLANST  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric tridiagonal matrix A.
!
!  Description
!  ===========
!
!  SLANST returns the value
!
!     SLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
!
!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
!
!  Arguments
!  =========
!
!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in SLANST as described
!          above.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, SLANST is
!          set to zero.
!
!  D       (input) REAL array, dimension (N)
!          The diagonal elements of A.
!
!  E       (input) REAL array, dimension (N-1)
!          The (n-1) sub-diagonal or super-diagonal elements of A.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      REAL               ANORM, SCALE, SUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLASSQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      IF( N<=0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
!
!        Find max(abs(A(i,j))).
!
         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            ANORM = MAX( ANORM, ABS( D( I ) ) )
            ANORM = MAX( ANORM, ABS( E( I ) ) )
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM=='1' .OR. &
               LSAME( NORM, 'I' ) ) THEN
!
!        Find norm1(A).
!
         IF( N==1 ) THEN
            ANORM = ABS( D( 1 ) )
         ELSE
            ANORM = MAX( ABS( D( 1 ) )+ABS( E( 1 ) ), &
                    ABS( E( N-1 ) )+ABS( D( N ) ) )
            DO 20 I = 2, N - 1
               ANORM = MAX( ANORM, ABS( D( I ) )+ABS( E( I ) )+ &
                       ABS( E( I-1 ) ) )
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
!
!        Find normF(A).
!
         SCALE = ZERO
         SUM = ONE
         IF( N>1 ) THEN
            CALL SLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         CALL SLASSQ( N, D, 1, SCALE, SUM )
         ANORM = SCALE*SQRT( SUM )
      END IF
!
      SLANST = ANORM
      RETURN
!
!     End of SLANST
!
      END
      SUBROUTINE SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      REAL               A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
!     ..
!
!  Purpose
!  =======
!
!  SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
!  matrix in standard form:
!
!       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
!       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
!
!  where either
!  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
!  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
!  conjugate eigenvalues.
!
!  Arguments
!  =========
!
!  A       (input/output) REAL
!  B       (input/output) REAL
!  C       (input/output) REAL
!  D       (input/output) REAL
!          On entry, the elements of the input matrix.
!          On exit, they are overwritten by the elements of the
!          standardised Schur form.
!
!  RT1R    (output) REAL
!  RT1I    (output) REAL
!  RT2R    (output) REAL
!  RT2I    (output) REAL
!          The real and imaginary parts of the eigenvalues. If the
!          eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the
!          eigenvalues are a complex conjugate pair, RT1I > 0.
!
!  CS      (output) REAL
!  SN      (output) REAL
!          Parameters of the rotation matrix.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0E+0, HALF = 0.5E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      REAL               AA, BB, CC, CS1, DD, P, SAB, SAC, SIGMA, SN1, &
                         TAU, TEMP
!     ..
!     .. External Functions ..
      REAL               SLAPY2
      EXTERNAL           SLAPY2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Initialize CS and SN
!
      CS = ONE
      SN = ZERO
!
      IF( C==ZERO ) THEN
         GO TO 10
!
      ELSE IF( B==ZERO ) THEN
!
!        Swap rows and columns
!
         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 10
      ELSE IF( (A-D)==ZERO .AND. SIGN( ONE, B )/= &
         SIGN( ONE, C ) ) THEN
         GO TO 10
      ELSE
!
!        Make diagonal elements equal
!
         TEMP = A - D
         P = HALF*TEMP
         SIGMA = B + C
         TAU = SLAPY2( SIGMA, TEMP )
         CS1 = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
         SN1 = -( P / ( TAU*CS1 ) )*SIGN( ONE, SIGMA )
!
!        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]
!                [ CC  DD ]   [ C  D ] [ SN1  CS1 ]
!
         AA = A*CS1 + B*SN1
         BB = -A*SN1 + B*CS1
         CC = C*CS1 + D*SN1
         DD = -C*SN1 + D*CS1
!
!        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]
!                [ C  D ]   [-SN1  CS1 ] [ CC  DD ]
!
         A = AA*CS1 + CC*SN1
         B = BB*CS1 + DD*SN1
         C = -AA*SN1 + CC*CS1
         D = -BB*SN1 + DD*CS1
!
!        Accumulate transformation
!
         TEMP = CS*CS1 - SN*SN1
         SN = CS*SN1 + SN*CS1
         CS = TEMP
!
         TEMP = HALF*( A+D )
         A = TEMP
         D = TEMP
!
         IF( C/=ZERO ) THEN
            IF ( B/=ZERO ) THEN
               IF( SIGN( ONE, B )==SIGN( ONE, C ) ) THEN
!
!                 Real eigenvalues: reduce to upper triangular form
!
                  SAB = SQRT( ABS( B ) )
                  SAC = SQRT( ABS( C ) )
                  P = SIGN( SAB*SAC, C )
                  TAU = ONE / SQRT( ABS( B+C ) )
                  A = TEMP + P
                  D = TEMP - P
                  B = B - C
                  C = ZERO
                  CS1 = SAB*TAU
                  SN1 = SAC*TAU
                  TEMP = CS*CS1 - SN*SN1
                  SN = CS*SN1 + SN*CS1
                  CS = TEMP
               END IF
            ELSE
               B = -C
               C = ZERO
               TEMP = CS
               CS = -SN
               SN = TEMP
            ENDIF
         ENDIF
      END IF
!
   10 CONTINUE
!
!     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
!
      RT1R = A
      RT2R = D
      IF( C==ZERO ) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
         RT2I = -RT1I
      END IF
      RETURN
!
!     End of SLANV2
!
      END
      SUBROUTINE SLAPTM( N, NRHS, ALPHA, D, E, X, LDX, BETA, B, LDB )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            LDB, LDX, N, NRHS
      REAL               ALPHA, BETA
!     ..
!     .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), E( * ), X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  SLAPTM multiplies an N by NRHS matrix X by a symmetric tridiagonal
!  matrix A and stores the result in a matrix B.  The operation has the
!  form
!
!     B := alpha * A * X + beta * B
!
!  where alpha may be either 1. or -1. and beta may be 0., 1., or -1.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrices X and B.
!
!  ALPHA   (input) REAL
!          The scalar alpha.  ALPHA must be 1. or -1.; otherwise,
!          it is assumed to be 0.
!
!  D       (input) REAL array, dimension (N)
!          The n diagonal elements of the tridiagonal matrix A.
!
!  E       (input) REAL array, dimension (N-1)
!          The (n-1) subdiagonal or superdiagonal elements of A.
!
!  X       (input) REAL array, dimension (LDX,NRHS)
!          The N by NRHS matrix X.
!
!  LDX     (input) INTEGER
!          The leading dimension of the array X.  LDX >= max(N,1).
!
!  BETA    (input) REAL
!          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,
!          it is assumed to be 1.
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the N by NRHS matrix B.
!          On exit, B is overwritten by the matrix expression
!          B := alpha * A * X + beta * B.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(N,1).
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. Executable Statements ..
!
      IF( N==0 ) &
         RETURN
!
!     Multiply B by BETA if BETA/=1.
!
      IF( BETA==ZERO ) THEN
         DO 20 J = 1, NRHS
            DO 10 I = 1, N
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE IF( BETA==-ONE ) THEN
         DO 40 J = 1, NRHS
            DO 30 I = 1, N
               B( I, J ) = -B( I, J )
   30       CONTINUE
   40    CONTINUE
      END IF
!
      IF( ALPHA==ONE ) THEN
!
!        Compute B := B + A*X
!
         DO 60 J = 1, NRHS
            IF( N==1 ) THEN
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + &
                           E( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) + E( N-1 )*X( N-1, J ) + &
                           D( N )*X( N, J )
               DO 50 I = 2, N - 1
                  B( I, J ) = B( I, J ) + E( I-1 )*X( I-1, J ) + &
                              D( I )*X( I, J ) + E( I )*X( I+1, J )
   50          CONTINUE
            END IF
   60    CONTINUE
      ELSE IF( ALPHA==-ONE ) THEN
!
!        Compute B := B - A*X
!
         DO 80 J = 1, NRHS
            IF( N==1 ) THEN
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
            ELSE
               B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - &
                           E( 1 )*X( 2, J )
               B( N, J ) = B( N, J ) - E( N-1 )*X( N-1, J ) - &
                           D( N )*X( N, J )
               DO 70 I = 2, N - 1
                  B( I, J ) = B( I, J ) - E( I-1 )*X( I-1, J ) - &
                              D( I )*X( I, J ) - E( I )*X( I+1, J )
   70          CONTINUE
            END IF
   80    CONTINUE
      END IF
      RETURN
!
!     End of SLAPTM
!
      END
      REAL             FUNCTION SLAPY2( X, Y )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               X, Y
!     ..
!
!  Purpose
!  =======
!
!  SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.
!
!  Arguments
!  =========
!
!  X       (input) REAL
!  Y       (input) REAL
!          X and Y specify the values x and y.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      REAL               W, XABS, YABS, Z
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      W = MAX( XABS, YABS )
      Z = MIN( XABS, YABS )
      IF( Z==ZERO ) THEN
         SLAPY2 = W
      ELSE
         SLAPY2 = W*SQRT( ONE+( Z / W )**2 )
      END IF
      RETURN
!
!     End of SLAPY2
!
      END
      REAL             FUNCTION SLAPY3( X, Y, Z )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      REAL               X, Y, Z
!     ..
!
!  Purpose
!  =======
!
!  SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!  unnecessary overflow.
!
!  Arguments
!  =========
!
!  X       (input) REAL
!  Y       (input) REAL
!  Z       (input) REAL
!          X, Y and Z specify the values x, y and z.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
!     ..
!     .. Local Scalars ..
      REAL               W, XABS, YABS, ZABS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W==ZERO ) THEN
         SLAPY3 = ZERO
      ELSE
         SLAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+ &
                  ( ZABS / W )**2 )
      END IF
      RETURN
!
!     End of SLAPY3
!
      END
      REAL FUNCTION SLARAN( ISEED )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
!     ..
!
!  Purpose
!  =======
!
!  SLARAN returns a random real number from a uniform (0,1)
!  distribution.
!
!  Arguments
!  =========
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  Further Details
!  ===============
!
!  This routine uses a multiplicative congruential method with modulus
!  2**48 and multiplier 33952834046453 (see G.S.Fishman,
!  'Multiplicative congruential random number generators with modulus
!  2**b: an exhaustive analysis for b = 32 and a partial analysis for
!  b = 48', Math. Comp. 189, pp 331-344, 1990).
!
!  48-bit integers are stored in 4 integer array elements with 12 bits
!  per element. Hence the routine is portable across machines with
!  integers of 32 bits or more.
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            M1, M2, M3, M4
      PARAMETER          ( M1 = 494, M2 = 322, M3 = 2508, M4 = 2549 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      INTEGER            IPW2
      REAL               R
      PARAMETER          ( IPW2 = 4096, R = ONE / IPW2 )
!     ..
!     .. Local Scalars ..
      INTEGER            IT1, IT2, IT3, IT4
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MOD, REAL
!     ..
!     .. Executable Statements ..
!
!     multiply the seed by the multiplier modulo 2**48
!
      IT4 = ISEED( 4 )*M4
      IT3 = IT4 / IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + ISEED( 3 )*M4 + ISEED( 4 )*M3
      IT2 = IT3 / IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + ISEED( 2 )*M4 + ISEED( 3 )*M3 + ISEED( 4 )*M2
      IT1 = IT2 / IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + ISEED( 1 )*M4 + ISEED( 2 )*M3 + ISEED( 3 )*M2 + &
            ISEED( 4 )*M1
      IT1 = MOD( IT1, IPW2 )
!
!     return updated seed
!
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
!
!     convert 48-bit integer to a real number in the interval (0,1)
!
      SLARAN = R*( REAL( IT1 )+R*( REAL( IT2 )+R*( REAL( IT3 )+R* &
               ( REAL( IT4 ) ) ) ) )
      RETURN
!
!     End of SLARAN
!
      END
      SUBROUTINE SLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      REAL               TAU
!     ..
!     .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) REAL array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.
!
!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.
!
!  TAU     (input) REAL
!          The value tau in the representation of H.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension
!                         (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C
!
         IF( TAU/=ZERO ) THEN
!
!           w := C' * v
!
            CALL SGEMV( 'Transpose', M, N, ONE, C, LDC, V, INCV, ZERO, &
                        WORK, 1 )
!
!           C := C - v * w'
!
            CALL SGER( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
!
!        Form  C * H
!
         IF( TAU/=ZERO ) THEN
!
!           w := C * v
!
            CALL SGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV, &
                        ZERO, WORK, 1 )
!
!           C := C - w * v'
!
            CALL SGER( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
!
!     End of SLARF
!
      END
      SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               ALPHA, TAU
!     ..
!     .. Array Arguments ..
      REAL               X( * )
!     ..
!
!  Purpose
!  =======
!
!  SLARFG generates a real elementary reflector H of order n, such
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
!  =========
!
!  N       (input) INTEGER
!          The order of the elementary reflector.
!
!  ALPHA   (input/output) REAL
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.
!
!  X       (input/output) REAL array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.
!
!  INCX    (input) INTEGER
!          The increment between elements of X. INCX > 0.
!
!  TAU     (output) REAL
!          The value tau.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J, KNT
      REAL               BETA, RSAFMN, SAFMIN, XNORM
!     ..
!     .. External Functions ..
      REAL               SLAMCH, SLAPY2, SNRM2
      EXTERNAL           SLAMCH, SLAPY2, SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL           SSCAL
!     ..
!     .. Executable Statements ..
!
      IF( N<=1 ) THEN
         TAU = ZERO
         RETURN
      END IF
!
      XNORM = SNRM2( N-1, X, INCX )
!
      IF( XNORM==ZERO ) THEN
!
!        H  =  I
!
         TAU = ZERO
      ELSE
!
!        general case
!
         BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
         SAFMIN = SLAMCH( 'S' ) / SLAMCH( 'E' )
         IF( ABS( BETA )<SAFMIN ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            RSAFMN = ONE / SAFMIN
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL SSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHA = ALPHA*RSAFMN
            IF( ABS( BETA )<SAFMIN ) &
               GO TO 10
!
!           New BETA is at most 1, at least SAFMIN
!
            XNORM = SNRM2( N-1, X, INCX )
            BETA = -SIGN( SLAPY2( ALPHA, XNORM ), ALPHA )
            TAU = ( BETA-ALPHA ) / BETA
            CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
!
!           If ALPHA is subnormal, it may lose relative accuracy
!
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = ( BETA-ALPHA ) / BETA
            CALL SSCAL( N-1, ONE / ( ALPHA-BETA ), X, INCX )
            ALPHA = BETA
         END IF
      END IF
!
      RETURN
!
!     End of SLARFG
!
      END
      SUBROUTINE SLARFX( SIDE, M, N, V, TAU, C, LDC, WORK )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            LDC, M, N
      REAL               TAU
!     ..
!     .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SLARFX applies a real elementary reflector H to a real m by n
!  matrix C, from either the left or the right. H is represented in the
!  form
!
!        H = I - tau * v * v'
!
!  where tau is a real scalar and v is a real vector.
!
!  If tau = 0, then H is taken to be the unit matrix
!
!  This version uses inline code if H has order < 11.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H
!
!  M       (input) INTEGER
!          The number of rows of the matrix C.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C.
!
!  V       (input) REAL array, dimension (M) if SIDE = 'L'
!                                     or (N) if SIDE = 'R'
!          The vector v in the representation of H.
!
!  TAU     (input) REAL
!          The value tau in the representation of H.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= (1,M).
!
!  WORK    (workspace) REAL array, dimension
!                      (N) if SIDE = 'L'
!                      or (M) if SIDE = 'R'
!          WORK is not referenced if H has order < 11.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            J
      REAL               SUM, T1, T10, T2, T3, T4, T5, T6, T7, T8, T9, &
                         V1, V10, V2, V3, V4, V5, V6, V7, V8, V9
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SGEMV, SGER
!     ..
!     .. Executable Statements ..
!
      IF( TAU==ZERO ) &
         RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  H * C, where H has order m.
!
         GO TO ( 10, 30, 50, 70, 90, 110, 130, 150, &
                 170, 190 )M
!
!        Code for general M
!
!        w := C'*v
!
         CALL SGEMV( 'Transpose', M, N, ONE, C, LDC, V, 1, ZERO, WORK, &
                     1 )
!
!        C := C - tau * v * w'
!
         CALL SGER( M, N, -TAU, V, 1, WORK, 1, C, LDC )
         GO TO 410
   10    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 20 J = 1, N
            C( 1, J ) = T1*C( 1, J )
   20    CONTINUE
         GO TO 410
   30    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 40 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
   40    CONTINUE
         GO TO 410
   50    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 60 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
   60    CONTINUE
         GO TO 410
   70    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 80 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
   80    CONTINUE
         GO TO 410
   90    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 100 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
  100    CONTINUE
         GO TO 410
  110    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 120 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
  120    CONTINUE
         GO TO 410
  130    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 140 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
  140    CONTINUE
         GO TO 410
  150    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 160 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
  160    CONTINUE
         GO TO 410
  170    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 180 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
  180    CONTINUE
         GO TO 410
  190    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 200 J = 1, N
            SUM = V1*C( 1, J ) + V2*C( 2, J ) + V3*C( 3, J ) + &
                  V4*C( 4, J ) + V5*C( 5, J ) + V6*C( 6, J ) + &
                  V7*C( 7, J ) + V8*C( 8, J ) + V9*C( 9, J ) + &
                  V10*C( 10, J )
            C( 1, J ) = C( 1, J ) - SUM*T1
            C( 2, J ) = C( 2, J ) - SUM*T2
            C( 3, J ) = C( 3, J ) - SUM*T3
            C( 4, J ) = C( 4, J ) - SUM*T4
            C( 5, J ) = C( 5, J ) - SUM*T5
            C( 6, J ) = C( 6, J ) - SUM*T6
            C( 7, J ) = C( 7, J ) - SUM*T7
            C( 8, J ) = C( 8, J ) - SUM*T8
            C( 9, J ) = C( 9, J ) - SUM*T9
            C( 10, J ) = C( 10, J ) - SUM*T10
  200    CONTINUE
         GO TO 410
      ELSE
!
!        Form  C * H, where H has order n.
!
         GO TO ( 210, 230, 250, 270, 290, 310, 330, 350, &
                 370, 390 )N
!
!        Code for general N
!
!        w := C * v
!
         CALL SGEMV( 'No transpose', M, N, ONE, C, LDC, V, 1, ZERO, &
                     WORK, 1 )
!
!        C := C - tau * w * v'
!
         CALL SGER( M, N, -TAU, WORK, 1, V, 1, C, LDC )
         GO TO 410
  210    CONTINUE
!
!        Special code for 1 x 1 Householder
!
         T1 = ONE - TAU*V( 1 )*V( 1 )
         DO 220 J = 1, M
            C( J, 1 ) = T1*C( J, 1 )
  220    CONTINUE
         GO TO 410
  230    CONTINUE
!
!        Special code for 2 x 2 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         DO 240 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
  240    CONTINUE
         GO TO 410
  250    CONTINUE
!
!        Special code for 3 x 3 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         DO 260 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
  260    CONTINUE
         GO TO 410
  270    CONTINUE
!
!        Special code for 4 x 4 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         DO 280 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
  280    CONTINUE
         GO TO 410
  290    CONTINUE
!
!        Special code for 5 x 5 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         DO 300 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
  300    CONTINUE
         GO TO 410
  310    CONTINUE
!
!        Special code for 6 x 6 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         DO 320 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
  320    CONTINUE
         GO TO 410
  330    CONTINUE
!
!        Special code for 7 x 7 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         DO 340 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
  340    CONTINUE
         GO TO 410
  350    CONTINUE
!
!        Special code for 8 x 8 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         DO 360 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
  360    CONTINUE
         GO TO 410
  370    CONTINUE
!
!        Special code for 9 x 9 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         DO 380 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
  380    CONTINUE
         GO TO 410
  390    CONTINUE
!
!        Special code for 10 x 10 Householder
!
         V1 = V( 1 )
         T1 = TAU*V1
         V2 = V( 2 )
         T2 = TAU*V2
         V3 = V( 3 )
         T3 = TAU*V3
         V4 = V( 4 )
         T4 = TAU*V4
         V5 = V( 5 )
         T5 = TAU*V5
         V6 = V( 6 )
         T6 = TAU*V6
         V7 = V( 7 )
         T7 = TAU*V7
         V8 = V( 8 )
         T8 = TAU*V8
         V9 = V( 9 )
         T9 = TAU*V9
         V10 = V( 10 )
         T10 = TAU*V10
         DO 400 J = 1, M
            SUM = V1*C( J, 1 ) + V2*C( J, 2 ) + V3*C( J, 3 ) + &
                  V4*C( J, 4 ) + V5*C( J, 5 ) + V6*C( J, 6 ) + &
                  V7*C( J, 7 ) + V8*C( J, 8 ) + V9*C( J, 9 ) + &
                  V10*C( J, 10 )
            C( J, 1 ) = C( J, 1 ) - SUM*T1
            C( J, 2 ) = C( J, 2 ) - SUM*T2
            C( J, 3 ) = C( J, 3 ) - SUM*T3
            C( J, 4 ) = C( J, 4 ) - SUM*T4
            C( J, 5 ) = C( J, 5 ) - SUM*T5
            C( J, 6 ) = C( J, 6 ) - SUM*T6
            C( J, 7 ) = C( J, 7 ) - SUM*T7
            C( J, 8 ) = C( J, 8 ) - SUM*T8
            C( J, 9 ) = C( J, 9 ) - SUM*T9
            C( J, 10 ) = C( J, 10 ) - SUM*T10
  400    CONTINUE
         GO TO 410
      END IF
  410 RETURN
!
!     End of SLARFX
!
      END
      REAL             FUNCTION SLARND( IDIST, ISEED )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IDIST
!     ..
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
!     ..
!
!  Purpose
!  =======
!
!  SLARND returns a random real number from a uniform or normal
!  distribution.
!
!  Arguments
!  =========
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  Further Details
!  ===============
!
!  This routine calls the auxiliary routine SLARAN to generate a random
!  real number from a uniform (0,1) distribution. The Box-Muller method
!  is used to transform numbers from a uniform to a normal distribution.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E+0, TWO = 2.0E+0 )
      REAL               TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663E+0 )
!     ..
!     .. Local Scalars ..
      REAL               T1, T2
!     ..
!     .. External Functions ..
      REAL               SLARAN
      EXTERNAL           SLARAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, SQRT
!     ..
!     .. Executable Statements ..
!
!     Generate a real random number from a uniform (0,1) distribution
!
      T1 = SLARAN( ISEED )
!
      IF( IDIST==1 ) THEN
!
!        uniform (0,1)
!
         SLARND = T1
      ELSE IF( IDIST==2 ) THEN
!
!        uniform (-1,1)
!
         SLARND = TWO*T1 - ONE
      ELSE IF( IDIST==3 ) THEN
!
!        normal (0,1)
!
         T2 = SLARAN( ISEED )
         SLARND = SQRT( -TWO*LOG( T1 ) )*COS( TWOPI*T2 )
      END IF
      RETURN
!
!     End of SLARND
!
      END
      SUBROUTINE SLARNV( IDIST, ISEED, N, X )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      INTEGER            IDIST, N
!     ..
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL               X( * )
!     ..
!
!  Purpose
!  =======
!
!  SLARNV returns a vector of n random real numbers from a uniform or
!  normal distribution.
!
!  Arguments
!  =========
!
!  IDIST   (input) INTEGER
!          Specifies the distribution of the random numbers:
!          = 1:  uniform (0,1)
!          = 2:  uniform (-1,1)
!          = 3:  normal (0,1)
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated.
!
!  X       (output) REAL array, dimension (N)
!          The generated random numbers.
!
!  Further Details
!  ===============
!
!  This routine calls the auxiliary routine SLARUV to generate random
!  real numbers from a uniform (0,1) distribution, in batches of up to
!  128 using vectorisable code. The Box-Muller method is used to
!  transform numbers from a uniform to a normal distribution.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, TWO
      PARAMETER          ( ONE = 1.0E+0, TWO = 2.0E+0 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      REAL               TWOPI
      PARAMETER          ( TWOPI = 6.2831853071795864769252867663E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, IL, IL2, IV
!     ..
!     .. Local Arrays ..
      REAL               U( LV )
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, MIN, SQRT
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARUV
!     ..
!     .. Executable Statements ..
!
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST==3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
!
!        Call SLARUV to generate IL2 numbers from a uniform (0,1)
!        distribution (IL2 <= LV)
!
         CALL SLARUV( ISEED, IL2, U )
!
         IF( IDIST==1 ) THEN
!
!           Copy generated numbers
!
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST==2 ) THEN
!
!           Convert generated numbers to uniform (-1,1) distribution
!
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST==3 ) THEN
!
!           Convert generated numbers to normal (0,1) distribution
!
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )* &
                             COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
!
!     End of SLARNV
!
      END
      SUBROUTINE SLARTG( F, G, CS, SN, R )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      REAL               CS, F, G, R, SN
!     ..
!
!  Purpose
!  =======
!
!  SLARTG generate a plane rotation so that
!
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]
!
!  This is a slower, more accurate version of the BLAS1 routine SROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in SBDSQR when
!        there are zeros on the diagonal).
!
!  If F exceeds G in magnitude, CS will be positive.
!
!  Arguments
!  =========
!
!  F       (input) REAL
!          The first component of vector to be rotated.
!
!  G       (input) REAL
!          The second component of vector to be rotated.
!
!  CS      (output) REAL
!          The cosine of the rotation.
!
!  SN      (output) REAL
!          The sine of the rotation.
!
!  R       (output) REAL
!          The nonzero component of the rotated vector.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      REAL               TWO
      PARAMETER          ( TWO = 2.0E0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      REAL               EPS, F1, G1, SAFMIN, SAFMN2, SAFMX2, SCALE
!     ..
!     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, INT, LOG, MAX, SQRT
!     ..
!     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
!     ..
!     .. Data statements ..
      DATA               FIRST / .TRUE. /
!     ..
!     .. Executable Statements ..
!
      IF( FIRST ) THEN
         FIRST = .FALSE.
         SAFMIN = SLAMCH( 'S' )
         EPS = SLAMCH( 'E' )
         SAFMN2 = SLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) / &
                  LOG( SLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      IF( G==ZERO ) THEN
         CS = ONE
         SN = ZERO
         R = F
      ELSE IF( F==ZERO ) THEN
         CS = ZERO
         SN = ONE
         R = G
      ELSE
         F1 = F
         G1 = G
         SCALE = MAX( ABS( F1 ), ABS( G1 ) )
         IF( SCALE>=SAFMX2 ) THEN
            COUNT = 0
   10       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE>=SAFMX2 ) &
               GO TO 10
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 20 I = 1, COUNT
               R = R*SAFMX2
   20       CONTINUE
         ELSE IF( SCALE<=SAFMN2 ) THEN
            COUNT = 0
   30       CONTINUE
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX( ABS( F1 ), ABS( G1 ) )
            IF( SCALE<=SAFMN2 ) &
               GO TO 30
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
            DO 40 I = 1, COUNT
               R = R*SAFMN2
   40       CONTINUE
         ELSE
            R = SQRT( F1**2+G1**2 )
            CS = F1 / R
            SN = G1 / R
         END IF
         IF( ABS( F )>ABS( G ) .AND. CS<ZERO ) THEN
            CS = -CS
            SN = -SN
            R = -R
         END IF
      END IF
      RETURN
!
!     End of SLARTG
!
      END
      SUBROUTINE SLARUV( ISEED, N, X )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            N
!     ..
!     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL               X( N )
!     ..
!
!  Purpose
!  =======
!
!  SLARUV returns a vector of n random real numbers from a uniform (0,1)
!  distribution (n <= 128).
!
!  This is an auxiliary routine called by SLARNV and CLARNV.
!
!  Arguments
!  =========
!
!  ISEED   (input/output) INTEGER array, dimension (4)
!          On entry, the seed of the random number generator; the array
!          elements must be between 0 and 4095, and ISEED(4) must be
!          odd.
!          On exit, the seed is updated.
!
!  N       (input) INTEGER
!          The number of random numbers to be generated. N <= 128.
!
!  X       (output) REAL array, dimension (N)
!          The generated random numbers.
!
!  Further Details
!  ===============
!
!  This routine uses a multiplicative congruential method with modulus
!  2**48 and multiplier 33952834046453 (see G.S.Fishman,
!  'Multiplicative congruential random number generators with modulus
!  2**b: an exhaustive analysis for b = 32 and a partial analysis for
!  b = 48', Math. Comp. 189, pp 331-344, 1990).
!
!  48-bit integers are stored in 4 integer array elements with 12 bits
!  per element. Hence the routine is portable across machines with
!  integers of 32 bits or more.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
      INTEGER            LV, IPW2
      REAL               R
      PARAMETER          ( LV = 128, IPW2 = 4096, R = ONE / IPW2 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, I1, I2, I3, I4, IT1, IT2, IT3, IT4, J
!     ..
!     .. Local Arrays ..
      INTEGER            MM( LV, 4 )
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN, MOD, REAL
!     ..
!     .. Data statements ..
      DATA               ( MM( 1, J ), J = 1, 4 ) / 494, 322, 2508, &
                         2549 /
      DATA               ( MM( 2, J ), J = 1, 4 ) / 2637, 789, 3754, &
                         1145 /
      DATA               ( MM( 3, J ), J = 1, 4 ) / 255, 1440, 1766, &
                         2253 /
      DATA               ( MM( 4, J ), J = 1, 4 ) / 2008, 752, 3572, &
                         305 /
      DATA               ( MM( 5, J ), J = 1, 4 ) / 1253, 2859, 2893, &
                         3301 /
      DATA               ( MM( 6, J ), J = 1, 4 ) / 3344, 123, 307, &
                         1065 /
      DATA               ( MM( 7, J ), J = 1, 4 ) / 4084, 1848, 1297, &
                         3133 /
      DATA               ( MM( 8, J ), J = 1, 4 ) / 1739, 643, 3966, &
                         2913 /
      DATA               ( MM( 9, J ), J = 1, 4 ) / 3143, 2405, 758, &
                         3285 /
      DATA               ( MM( 10, J ), J = 1, 4 ) / 3468, 2638, 2598, &
                         1241 /
      DATA               ( MM( 11, J ), J = 1, 4 ) / 688, 2344, 3406, &
                         1197 /
      DATA               ( MM( 12, J ), J = 1, 4 ) / 1657, 46, 2922, &
                         3729 /
      DATA               ( MM( 13, J ), J = 1, 4 ) / 1238, 3814, 1038, &
                         2501 /
      DATA               ( MM( 14, J ), J = 1, 4 ) / 3166, 913, 2934, &
                         1673 /
      DATA               ( MM( 15, J ), J = 1, 4 ) / 1292, 3649, 2091, &
                         541 /
      DATA               ( MM( 16, J ), J = 1, 4 ) / 3422, 339, 2451, &
                         2753 /
      DATA               ( MM( 17, J ), J = 1, 4 ) / 1270, 3808, 1580, &
                         949 /
      DATA               ( MM( 18, J ), J = 1, 4 ) / 2016, 822, 1958, &
                         2361 /
      DATA               ( MM( 19, J ), J = 1, 4 ) / 154, 2832, 2055, &
                         1165 /
      DATA               ( MM( 20, J ), J = 1, 4 ) / 2862, 3078, 1507, &
                         4081 /
      DATA               ( MM( 21, J ), J = 1, 4 ) / 697, 3633, 1078, &
                         2725 /
      DATA               ( MM( 22, J ), J = 1, 4 ) / 1706, 2970, 3273, &
                         3305 /
      DATA               ( MM( 23, J ), J = 1, 4 ) / 491, 637, 17, &
                         3069 /
      DATA               ( MM( 24, J ), J = 1, 4 ) / 931, 2249, 854, &
                         3617 /
      DATA               ( MM( 25, J ), J = 1, 4 ) / 1444, 2081, 2916, &
                         3733 /
      DATA               ( MM( 26, J ), J = 1, 4 ) / 444, 4019, 3971, &
                         409 /
      DATA               ( MM( 27, J ), J = 1, 4 ) / 3577, 1478, 2889, &
                         2157 /
      DATA               ( MM( 28, J ), J = 1, 4 ) / 3944, 242, 3831, &
                         1361 /
      DATA               ( MM( 29, J ), J = 1, 4 ) / 2184, 481, 2621, &
                         3973 /
      DATA               ( MM( 30, J ), J = 1, 4 ) / 1661, 2075, 1541, &
                         1865 /
      DATA               ( MM( 31, J ), J = 1, 4 ) / 3482, 4058, 893, &
                         2525 /
      DATA               ( MM( 32, J ), J = 1, 4 ) / 657, 622, 736, &
                         1409 /
      DATA               ( MM( 33, J ), J = 1, 4 ) / 3023, 3376, 3992, &
                         3445 /
      DATA               ( MM( 34, J ), J = 1, 4 ) / 3618, 812, 787, &
                         3577 /
      DATA               ( MM( 35, J ), J = 1, 4 ) / 1267, 234, 2125, &
                         77 /
      DATA               ( MM( 36, J ), J = 1, 4 ) / 1828, 641, 2364, &
                         3761 /
      DATA               ( MM( 37, J ), J = 1, 4 ) / 164, 4005, 2460, &
                         2149 /
      DATA               ( MM( 38, J ), J = 1, 4 ) / 3798, 1122, 257, &
                         1449 /
      DATA               ( MM( 39, J ), J = 1, 4 ) / 3087, 3135, 1574, &
                         3005 /
      DATA               ( MM( 40, J ), J = 1, 4 ) / 2400, 2640, 3912, &
                         225 /
      DATA               ( MM( 41, J ), J = 1, 4 ) / 2870, 2302, 1216, &
                         85 /
      DATA               ( MM( 42, J ), J = 1, 4 ) / 3876, 40, 3248, &
                         3673 /
      DATA               ( MM( 43, J ), J = 1, 4 ) / 1905, 1832, 3401, &
                         3117 /
      DATA               ( MM( 44, J ), J = 1, 4 ) / 1593, 2247, 2124, &
                         3089 /
      DATA               ( MM( 45, J ), J = 1, 4 ) / 1797, 2034, 2762, &
                         1349 /
      DATA               ( MM( 46, J ), J = 1, 4 ) / 1234, 2637, 149, &
                         2057 /
      DATA               ( MM( 47, J ), J = 1, 4 ) / 3460, 1287, 2245, &
                         413 /
      DATA               ( MM( 48, J ), J = 1, 4 ) / 328, 1691, 166, &
                         65 /
      DATA               ( MM( 49, J ), J = 1, 4 ) / 2861, 496, 466, &
                         1845 /
      DATA               ( MM( 50, J ), J = 1, 4 ) / 1950, 1597, 4018, &
                         697 /
      DATA               ( MM( 51, J ), J = 1, 4 ) / 617, 2394, 1399, &
                         3085 /
      DATA               ( MM( 52, J ), J = 1, 4 ) / 2070, 2584, 190, &
                         3441 /
      DATA               ( MM( 53, J ), J = 1, 4 ) / 3331, 1843, 2879, &
                         1573 /
      DATA               ( MM( 54, J ), J = 1, 4 ) / 769, 336, 153, &
                         3689 /
      DATA               ( MM( 55, J ), J = 1, 4 ) / 1558, 1472, 2320, &
                         2941 /
      DATA               ( MM( 56, J ), J = 1, 4 ) / 2412, 2407, 18, &
                         929 /
      DATA               ( MM( 57, J ), J = 1, 4 ) / 2800, 433, 712, &
                         533 /
      DATA               ( MM( 58, J ), J = 1, 4 ) / 189, 2096, 2159, &
                         2841 /
      DATA               ( MM( 59, J ), J = 1, 4 ) / 287, 1761, 2318, &
                         4077 /
      DATA               ( MM( 60, J ), J = 1, 4 ) / 2045, 2810, 2091, &
                         721 /
      DATA               ( MM( 61, J ), J = 1, 4 ) / 1227, 566, 3443, &
                         2821 /
      DATA               ( MM( 62, J ), J = 1, 4 ) / 2838, 442, 1510, &
                         2249 /
      DATA               ( MM( 63, J ), J = 1, 4 ) / 209, 41, 449, &
                         2397 /
      DATA               ( MM( 64, J ), J = 1, 4 ) / 2770, 1238, 1956, &
                         2817 /
      DATA               ( MM( 65, J ), J = 1, 4 ) / 3654, 1086, 2201, &
                         245 /
      DATA               ( MM( 66, J ), J = 1, 4 ) / 3993, 603, 3137, &
                         1913 /
      DATA               ( MM( 67, J ), J = 1, 4 ) / 192, 840, 3399, &
                         1997 /
      DATA               ( MM( 68, J ), J = 1, 4 ) / 2253, 3168, 1321, &
                         3121 /
      DATA               ( MM( 69, J ), J = 1, 4 ) / 3491, 1499, 2271, &
                         997 /
      DATA               ( MM( 70, J ), J = 1, 4 ) / 2889, 1084, 3667, &
                         1833 /
      DATA               ( MM( 71, J ), J = 1, 4 ) / 2857, 3438, 2703, &
                         2877 /
      DATA               ( MM( 72, J ), J = 1, 4 ) / 2094, 2408, 629, &
                         1633 /
      DATA               ( MM( 73, J ), J = 1, 4 ) / 1818, 1589, 2365, &
                         981 /
      DATA               ( MM( 74, J ), J = 1, 4 ) / 688, 2391, 2431, &
                         2009 /
      DATA               ( MM( 75, J ), J = 1, 4 ) / 1407, 288, 1113, &
                         941 /
      DATA               ( MM( 76, J ), J = 1, 4 ) / 634, 26, 3922, &
                         2449 /
      DATA               ( MM( 77, J ), J = 1, 4 ) / 3231, 512, 2554, &
                         197 /
      DATA               ( MM( 78, J ), J = 1, 4 ) / 815, 1456, 184, &
                         2441 /
      DATA               ( MM( 79, J ), J = 1, 4 ) / 3524, 171, 2099, &
                         285 /
      DATA               ( MM( 80, J ), J = 1, 4 ) / 1914, 1677, 3228, &
                         1473 /
      DATA               ( MM( 81, J ), J = 1, 4 ) / 516, 2657, 4012, &
                         2741 /
      DATA               ( MM( 82, J ), J = 1, 4 ) / 164, 2270, 1921, &
                         3129 /
      DATA               ( MM( 83, J ), J = 1, 4 ) / 303, 2587, 3452, &
                         909 /
      DATA               ( MM( 84, J ), J = 1, 4 ) / 2144, 2961, 3901, &
                         2801 /
      DATA               ( MM( 85, J ), J = 1, 4 ) / 3480, 1970, 572, &
                         421 /
      DATA               ( MM( 86, J ), J = 1, 4 ) / 119, 1817, 3309, &
                         4073 /
      DATA               ( MM( 87, J ), J = 1, 4 ) / 3357, 676, 3171, &
                         2813 /
      DATA               ( MM( 88, J ), J = 1, 4 ) / 837, 1410, 817, &
                         2337 /
      DATA               ( MM( 89, J ), J = 1, 4 ) / 2826, 3723, 3039, &
                         1429 /
      DATA               ( MM( 90, J ), J = 1, 4 ) / 2332, 2803, 1696, &
                         1177 /
      DATA               ( MM( 91, J ), J = 1, 4 ) / 2089, 3185, 1256, &
                         1901 /
      DATA               ( MM( 92, J ), J = 1, 4 ) / 3780, 184, 3715, &
                         81 /
      DATA               ( MM( 93, J ), J = 1, 4 ) / 1700, 663, 2077, &
                         1669 /
      DATA               ( MM( 94, J ), J = 1, 4 ) / 3712, 499, 3019, &
                         2633 /
      DATA               ( MM( 95, J ), J = 1, 4 ) / 150, 3784, 1497, &
                         2269 /
      DATA               ( MM( 96, J ), J = 1, 4 ) / 2000, 1631, 1101, &
                         129 /
      DATA               ( MM( 97, J ), J = 1, 4 ) / 3375, 1925, 717, &
                         1141 /
      DATA               ( MM( 98, J ), J = 1, 4 ) / 1621, 3912, 51, &
                         249 /
      DATA               ( MM( 99, J ), J = 1, 4 ) / 3090, 1398, 981, &
                         3917 /
      DATA               ( MM( 100, J ), J = 1, 4 ) / 3765, 1349, 1978, &
                         2481 /
      DATA               ( MM( 101, J ), J = 1, 4 ) / 1149, 1441, 1813, &
                         3941 /
      DATA               ( MM( 102, J ), J = 1, 4 ) / 3146, 2224, 3881, &
                         2217 /
      DATA               ( MM( 103, J ), J = 1, 4 ) / 33, 2411, 76, &
                         2749 /
      DATA               ( MM( 104, J ), J = 1, 4 ) / 3082, 1907, 3846, &
                         3041 /
      DATA               ( MM( 105, J ), J = 1, 4 ) / 2741, 3192, 3694, &
                         1877 /
      DATA               ( MM( 106, J ), J = 1, 4 ) / 359, 2786, 1682, &
                         345 /
      DATA               ( MM( 107, J ), J = 1, 4 ) / 3316, 382, 124, &
                         2861 /
      DATA               ( MM( 108, J ), J = 1, 4 ) / 1749, 37, 1660, &
                         1809 /
      DATA               ( MM( 109, J ), J = 1, 4 ) / 185, 759, 3997, &
                         3141 /
      DATA               ( MM( 110, J ), J = 1, 4 ) / 2784, 2948, 479, &
                         2825 /
      DATA               ( MM( 111, J ), J = 1, 4 ) / 2202, 1862, 1141, &
                         157 /
      DATA               ( MM( 112, J ), J = 1, 4 ) / 2199, 3802, 886, &
                         2881 /
      DATA               ( MM( 113, J ), J = 1, 4 ) / 1364, 2423, 3514, &
                         3637 /
      DATA               ( MM( 114, J ), J = 1, 4 ) / 1244, 2051, 1301, &
                         1465 /
      DATA               ( MM( 115, J ), J = 1, 4 ) / 2020, 2295, 3604, &
                         2829 /
      DATA               ( MM( 116, J ), J = 1, 4 ) / 3160, 1332, 1888, &
                         2161 /
      DATA               ( MM( 117, J ), J = 1, 4 ) / 2785, 1832, 1836, &
                         3365 /
      DATA               ( MM( 118, J ), J = 1, 4 ) / 2772, 2405, 1990, &
                         361 /
      DATA               ( MM( 119, J ), J = 1, 4 ) / 1217, 3638, 2058, &
                         2685 /
      DATA               ( MM( 120, J ), J = 1, 4 ) / 1822, 3661, 692, &
                         3745 /
      DATA               ( MM( 121, J ), J = 1, 4 ) / 1245, 327, 1194, &
                         2325 /
      DATA               ( MM( 122, J ), J = 1, 4 ) / 2252, 3660, 20, &
                         3609 /
      DATA               ( MM( 123, J ), J = 1, 4 ) / 3904, 716, 3285, &
                         3821 /
      DATA               ( MM( 124, J ), J = 1, 4 ) / 2774, 1842, 2046, &
                         3537 /
      DATA               ( MM( 125, J ), J = 1, 4 ) / 997, 3987, 2107, &
                         517 /
      DATA               ( MM( 126, J ), J = 1, 4 ) / 2573, 1368, 3508, &
                         3017 /
      DATA               ( MM( 127, J ), J = 1, 4 ) / 1148, 1848, 3525, &
                         2141 /
      DATA               ( MM( 128, J ), J = 1, 4 ) / 545, 2366, 3801, &
                         1537 /
!     ..
!     .. Executable Statements ..
!
      I1 = ISEED( 1 )
      I2 = ISEED( 2 )
      I3 = ISEED( 3 )
      I4 = ISEED( 4 )
!
      DO 10 I = 1, MIN( N, LV )
!
!        Multiply the seed by i-th power of the multiplier modulo 2**48
!
         IT4 = I4*MM( I, 4 )
         IT3 = IT4 / IPW2
         IT4 = IT4 - IPW2*IT3
         IT3 = IT3 + I3*MM( I, 4 ) + I4*MM( I, 3 )
         IT2 = IT3 / IPW2
         IT3 = IT3 - IPW2*IT2
         IT2 = IT2 + I2*MM( I, 4 ) + I3*MM( I, 3 ) + I4*MM( I, 2 )
         IT1 = IT2 / IPW2
         IT2 = IT2 - IPW2*IT1
         IT1 = IT1 + I1*MM( I, 4 ) + I2*MM( I, 3 ) + I3*MM( I, 2 ) + &
               I4*MM( I, 1 )
         IT1 = MOD( IT1, IPW2 )
!
!        Convert 48-bit integer to a real number in the interval (0,1)
!
         X( I ) = R*( REAL( IT1 )+R*( REAL( IT2 )+R*( REAL( IT3 )+R* &
                  REAL( IT4 ) ) ) )
   10 CONTINUE
!
!     Return final value of seed
!
      ISEED( 1 ) = IT1
      ISEED( 2 ) = IT2
      ISEED( 3 ) = IT3
      ISEED( 4 ) = IT4
      RETURN
!
!     End of SLARUV
!
      END
      SUBROUTINE SLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      REAL               CFROM, CTO
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASCL multiplies the M by N real matrix A by the real scalar
!  CTO/CFROM.  This is done without over/underflow as long as the final
!  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
!  A may be full, upper triangular, lower triangular, upper Hessenberg,
!  or banded.
!
!  Arguments
!  =========
!
!  TYPE    (input) CHARACTER*1
!          TYPE indices the storage type of the input matrix.
!          = 'G':  A is a full matrix.
!          = 'L':  A is a lower triangular matrix.
!          = 'U':  A is an upper triangular matrix.
!          = 'H':  A is an upper Hessenberg matrix.
!          = 'B':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the lower
!                  half stored.
!          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
!                  and upper bandwidth KU and with the only the upper
!                  half stored.
!          = 'Z':  A is a band matrix with lower bandwidth KL and upper
!                  bandwidth KU.
!
!  KL      (input) INTEGER
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  KU      (input) INTEGER
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',
!          'Q' or 'Z'.
!
!  CFROM   (input) REAL
!  CTO     (input) REAL
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
!          without over/underflow if the final result CTO*A(I,J)/CFROM
!          can be represented without over/underflow.  CFROM must be
!          nonzero.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) REAL array, dimension (LDA,M)
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
!          storage type.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  INFO    (output) INTEGER
!          0  - successful exit
!          <0 - if INFO = -i, the i-th argument had an illegal value.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH
      EXTERNAL           LSAME, SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
!
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
!
      IF( ITYPE==-1 ) THEN
         INFO = -1
      ELSE IF( CFROM==ZERO ) THEN
         INFO = -4
      ELSE IF( M<0 ) THEN
         INFO = -6
      ELSE IF( N<0 .OR. ( ITYPE==4 .AND. N/=M ) .OR. &
               ( ITYPE==5 .AND. N/=M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE<=3 .AND. LDA<MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE>=4 ) THEN
         IF( KL<0 .OR. KL>MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU<0 .OR. KU>MAX( N-1, 0 ) .OR. &
                  ( ( ITYPE==4 .OR. ITYPE==5 ) .AND. KL/=KU ) ) &
                   THEN
            INFO = -3
         ELSE IF( ( ITYPE==4 .AND. LDA<KL+1 ) .OR. &
                  ( ITYPE==5 .AND. LDA<KU+1 ) .OR. &
                  ( ITYPE==6 .AND. LDA<2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
!
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SLASCL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N==0 .OR. M==0 ) &
         RETURN
!
!     Get machine parameters
!
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 )>ABS( CTOC ) .AND. CTOC/=ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 )>ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
!
      IF( ITYPE==0 ) THEN
!
!        Full matrix
!
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
!
      ELSE IF( ITYPE==1 ) THEN
!
!        Lower triangular matrix
!
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
!
      ELSE IF( ITYPE==2 ) THEN
!
!        Upper triangular matrix
!
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
!
      ELSE IF( ITYPE==3 ) THEN
!
!        Upper Hessenberg matrix
!
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
!
      ELSE IF( ITYPE==4 ) THEN
!
!        Lower half of a symmetric band matrix
!
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
!
      ELSE IF( ITYPE==5 ) THEN
!
!        Upper half of a symmetric band matrix
!
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( ITYPE==6 ) THEN
!
!        Band matrix
!
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
!
      END IF
!
      IF( .NOT.DONE ) &
         GO TO 10
!
      RETURN
!
!     End of SLASCL
!
      END
      SUBROUTINE SLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      REAL               ALPHA, BETA
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  ALPHA   (input) REAL
!          The constant to which the offdiagonal elements are to be set.
!
!  BETA    (input) REAL
!          The constant to which the diagonal elements are to be set.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:
!
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MIN
!     ..
!     .. Executable Statements ..
!
      IF( LSAME( UPLO, 'U' ) ) THEN
!
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
!
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO 40 J = 1, MIN( M, N )
            DO 30 I = J + 1, M
               A( I, J ) = ALPHA
   30       CONTINUE
   40    CONTINUE
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO 60 J = 1, N
            DO 50 I = 1, M
               A( I, J ) = ALPHA
   50       CONTINUE
   60    CONTINUE
      END IF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO 70 I = 1, MIN( M, N )
         A( I, I ) = BETA
   70 CONTINUE
!
      RETURN
!
!     End of SLASET
!
      END
      SUBROUTINE SLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), C( * ), S( * )
!     ..
!
!  Purpose
!  =======
!
!  SLASR   performs the transformation
!
!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
!
!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
!
!  where A is an m by n real matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
!  and z = n when SIDE = 'R' or 'r' ):
!
!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
!
!     P = P( z - 1 )*...*P( 2 )*P( 1 ),
!
!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
!
!     P = P( 1 )*P( 2 )*...*P( z - 1 ),
!
!  where  P( k ) is a plane rotation matrix for the following planes:
!
!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
!        the plane ( k, k + 1 )
!
!     when  PIVOT = 'T' or 't'  ( Top pivot ),
!        the plane ( 1, k + 1 )
!
!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
!        the plane ( k, z )
!
!  c( k ) and s( k )  must contain the  cosine and sine that define the
!  matrix  P( k ).  The two by two plane rotation part of the matrix
!  P( k ), R( k ), is assumed to be of the form
!
!     R( k ) = (  c( k )  s( k ) ).
!              ( -s( k )  c( k ) )
!
!  This version vectorises across rows of the array A when SIDE = 'L'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'
!
!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
!
!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.
!
!  C, S    (input) REAL arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R( k ) = (  c( k )  s( k ) ).
!                   ( -s( k )  c( k ) )
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, INFO, J
      REAL               CTEMP, STEMP, TEMP
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT, &
               'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) ) &
                THEN
         INFO = 3
      ELSE IF( M<0 ) THEN
         INFO = 4
      ELSE IF( N<0 ) THEN
         INFO = 5
      ELSE IF( LDA<MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SLASR ', INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( ( M==0 ) .OR. ( N==0 ) ) &
         RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
!
!        Form  P * A
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
!
!        Form A * P'
!
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP/=ONE ) .OR. ( STEMP/=ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
!
      RETURN
!
!     End of SLASR
!
      END
      SUBROUTINE SLASRT( ID, N, D, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      REAL               D( * )
!     ..
!
!  Purpose
!  =======
!
!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).
!
!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.
!
!  Arguments
!  =========
!
!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.
!
!  N       (input) INTEGER
!          The length of the array D.
!
!  D       (input/output) REAL array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
!     ..
!     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT
      REAL               D1, D2, D3, DMNMX, TMP
!     ..
!     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input paramters.
!
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR==-1 ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SLASRT', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N<=1 ) &
         RETURN
!
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START<=SELECT .AND. ENDD-START>0 ) THEN
!
!        Do Insertion sort on D( START:ENDD )
!
         IF( DIR==0 ) THEN
!
!           Sort into decreasing order
!
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( D( J )>D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
!
         ELSE
!
!           Sort into increasing order
!
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( D( J )<D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
!
         END IF
!
      ELSE IF( ENDD-START>SELECT ) THEN
!
!        Partition D( START:ENDD ) and stack parts, largest one first
!
!        Choose partition entry as median of 3
!
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF( D1<D2 ) THEN
            IF( D3<D1 ) THEN
               DMNMX = D1
            ELSE IF( D3<D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            END IF
         ELSE
            IF( D3<D2 ) THEN
               DMNMX = D2
            ELSE IF( D3<D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            END IF
         END IF
!
         IF( DIR==0 ) THEN
!
!           Sort into decreasing order
!
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( D( J )<DMNMX ) &
               GO TO 70
   80       CONTINUE
            I = I + 1
            IF( D( I )>DMNMX ) &
               GO TO 80
            IF( I<J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            END IF
            IF( J-START>ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
!
!           Sort into increasing order
!
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( D( J )>DMNMX ) &
               GO TO 100
  110       CONTINUE
            I = I + 1
            IF( D( I )<DMNMX ) &
               GO TO 110
            IF( I<J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            END IF
            IF( J-START>ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT>0 ) &
         GO TO 10
      RETURN
!
!     End of SLASRT
!
      END
      SUBROUTINE SLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, N
      REAL               SCALE, SUMSQ
!     ..
!     .. Array Arguments ..
      REAL               X( * )
!     ..
!
!  Purpose
!  =======
!
!  SLASSQ  returns the values  scl  and  smsq  such that
!
!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!
!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value
!
!     scl = max( scale, abs( x( i ) ) ).
!
!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.
!
!  The routine makes only one pass through the vector x.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of elements to be used from the vector X.
!
!  X       (input) REAL
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!
!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.
!
!  SCALE   (input/output) REAL
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.
!
!  SUMSQ   (input/output) REAL
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            IX
      REAL               ABSXI
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS
!     ..
!     .. Executable Statements ..
!
      IF( N>0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( X( IX )/=ZERO ) THEN
               ABSXI = ABS( X( IX ) )
               IF( SCALE<ABSXI ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / ABSXI )**2
                  SCALE = ABSXI
               ELSE
                  SUMSQ = SUMSQ + ( ABSXI / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
      RETURN
!
!     End of SLASSQ
!
      END
      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
!     ..
!     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.
!
!  A       (input/output) REAL array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!
!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.
!
!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.
!
!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.
!
!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, IP, IX
!     ..
!     .. External Subroutines ..
      EXTERNAL           SSWAP
!     ..
!     .. Executable Statements ..
!
!     Interchange row I with row IPIV(I) for each of rows K1 through K2.
!
      IF( INCX==0 ) &
         RETURN
      IF( INCX>0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX==1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP/=I ) &
               CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX>1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP/=I ) &
               CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX<0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP/=I ) &
               CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
!
      RETURN
!
!     End of SLASWP
!
      END
      SUBROUTINE SLASY2( LTRANL, LTRANR, ISGN, N1, N2, TL, LDTL, TR, &
                         LDTR, B, LDB, SCALE, X, LDX, XNORM, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     October 31, 1992
!
!     .. Scalar Arguments ..
      LOGICAL            LTRANL, LTRANR
      INTEGER            INFO, ISGN, LDB, LDTL, LDTR, LDX, N1, N2
      REAL               SCALE, XNORM
!     ..
!     .. Array Arguments ..
      REAL               B( LDB, * ), TL( LDTL, * ), TR( LDTR, * ), &
                         X( LDX, * )
!     ..
!
!  Purpose
!  =======
!
!  SLASY2 solves for the N1 by N2 matrix X, 1 <= N1,N2 <= 2, in
!
!         op(TL)*X + ISGN*X*op(TR) = SCALE*B,
!
!  where TL is N1 by N1, TR is N2 by N2, B is N1 by N2, and ISGN = 1 or
!  -1.  op(T) = T or T', where T' denotes the transpose of T.
!
!  Arguments
!  =========
!
!  LTRANL  (input) LOGICAL
!          On entry, LTRANL specifies the op(TL):
!             = .FALSE., op(TL) = TL,
!             = .TRUE., op(TL) = TL'.
!
!  LTRANR  (input) LOGICAL
!          On entry, LTRANR specifies the op(TR):
!            = .FALSE., op(TR) = TR,
!            = .TRUE., op(TR) = TR'.
!
!  ISGN    (input) INTEGER
!          On entry, ISGN specifies the sign of the equation
!          as described before. ISGN may only be 1 or -1.
!
!  N1      (input) INTEGER
!          On entry, N1 specifies the order of matrix TL.
!          N1 may only be 0, 1 or 2.
!
!  N2      (input) INTEGER
!          On entry, N2 specifies the order of matrix TR.
!          N2 may only be 0, 1 or 2.
!
!  TL      (input) REAL array, dimension (LDTL,2)
!          On entry, TL contains an N1 by N1 matrix.
!
!  LDTL    (input) INTEGER
!          The leading dimension of the matrix TL. LDTL >= max(1,N1).
!
!  TR      (input) REAL array, dimension (LDTR,2)
!          On entry, TR contains an N2 by N2 matrix.
!
!  LDTR    (input) INTEGER
!          The leading dimension of the matrix TR. LDTR >= max(1,N2).
!
!  B       (input) REAL array, dimension (LDB,2)
!          On entry, the N1 by N2 matrix B contains the right-hand
!          side of the equation.
!
!  LDB     (input) INTEGER
!          The leading dimension of the matrix B. LDB >= max(1,N1).
!
!  SCALE   (output) REAL
!          On exit, SCALE contains the scale factor. SCALE is chosen
!          less than or equal to 1 to prevent the solution overflowing.
!
!  X       (output) REAL array, dimension (LDX,2)
!          On exit, X contains the N1 by N2 solution.
!
!  LDX     (input) INTEGER
!          The leading dimension of the matrix X. LDX >= max(1,N1).
!
!  XNORM   (output) REAL
!          On exit, XNORM is the infinity-norm of the solution.
!
!  INFO    (output) INTEGER
!          On exit, INFO is set to
!             0: successful exit.
!             1: TL and TR have too close eigenvalues, so TL or
!                TR is perturbed to get a nonsingular equation.
!          NOTE: In the interests of speed, this routine does not
!                check the inputs for errors.
!
! =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      REAL               TWO, HALF, EIGHT
      PARAMETER          ( TWO = 2.0E+0, HALF = 0.5E+0, EIGHT = 8.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            BSWAP, XSWAP
      INTEGER            I, IP, IPIV, IPSV, J, JP, JPSV, K
      REAL               BET, EPS, GAM, L21, SGN, SMIN, SMLNUM, TAU1, &
                         TEMP, U11, U12, U22, XMAX
!     ..
!     .. Local Arrays ..
      LOGICAL            BSWPIV( 4 ), XSWPIV( 4 )
      INTEGER            JPIV( 4 ), LOCL21( 4 ), LOCU12( 4 ), &
                         LOCU22( 4 )
      REAL               BTMP( 4 ), T16( 4, 4 ), TMP( 4 ), X2( 2 )
!     ..
!     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SLAMCH
      EXTERNAL           ISAMAX, SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           SCOPY, SSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
!     ..
!     .. Data statements ..
      DATA               LOCU12 / 3, 4, 1, 2 / , LOCL21 / 2, 1, 4, 3 / , &
                         LOCU22 / 4, 3, 2, 1 /
      DATA               XSWPIV / .FALSE., .FALSE., .TRUE., .TRUE. /
      DATA               BSWPIV / .FALSE., .TRUE., .FALSE., .TRUE. /
!     ..
!     .. Executable Statements ..
!
!     Do not check the input parameters for errors
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N1==0 .OR. N2==0 ) &
         RETURN
!
!     Set constants to control overflow
!
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      SGN = ISGN
!
      K = N1 + N1 + N2 - 2
      GO TO ( 10, 20, 30, 50 )K
!
!     1 by 1: TL11*X + SGN*X*TR11 = B11
!
   10 CONTINUE
      TAU1 = TL( 1, 1 ) + SGN*TR( 1, 1 )
      BET = ABS( TAU1 )
      IF( BET<=SMLNUM ) THEN
         TAU1 = SMLNUM
         BET = SMLNUM
         INFO = 1
      END IF
!
      SCALE = ONE
      GAM = ABS( B( 1, 1 ) )
      IF( SMLNUM*GAM>BET ) &
         SCALE = ONE / GAM
!
      X( 1, 1 ) = ( B( 1, 1 )*SCALE ) / TAU1
      XNORM = ABS( X( 1, 1 ) )
      RETURN
!
!     1 by 2:
!     TL11*[X11 X12] + ISGN*[X11 X12]*op[TR11 TR12]  = [B11 B12]
!                                       [TR21 TR22]
!
   20 CONTINUE
!
      SMIN = MAX( EPS*MAX( ABS( TL( 1, 1 ) ), ABS( TR( 1, 1 ) ), &
             ABS( TR( 1, 2 ) ), ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) ), &
             SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      IF( LTRANR ) THEN
         TMP( 2 ) = SGN*TR( 2, 1 )
         TMP( 3 ) = SGN*TR( 1, 2 )
      ELSE
         TMP( 2 ) = SGN*TR( 1, 2 )
         TMP( 3 ) = SGN*TR( 2, 1 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 1, 2 )
      GO TO 40
!
!     2 by 1:
!          op[TL11 TL12]*[X11] + ISGN* [X11]*TR11  = [B11]
!            [TL21 TL22] [X21]         [X21]         [B21]
!
   30 CONTINUE
      SMIN = MAX( EPS*MAX( ABS( TR( 1, 1 ) ), ABS( TL( 1, 1 ) ), &
             ABS( TL( 1, 2 ) ), ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) ), &
             SMLNUM )
      TMP( 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      TMP( 4 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      IF( LTRANL ) THEN
         TMP( 2 ) = TL( 1, 2 )
         TMP( 3 ) = TL( 2, 1 )
      ELSE
         TMP( 2 ) = TL( 2, 1 )
         TMP( 3 ) = TL( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
   40 CONTINUE
!
!     Solve 2 by 2 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
      IPIV = ISAMAX( 4, TMP, 1 )
      U11 = TMP( IPIV )
      IF( ABS( U11 )<=SMIN ) THEN
         INFO = 1
         U11 = SMIN
      END IF
      U12 = TMP( LOCU12( IPIV ) )
      L21 = TMP( LOCL21( IPIV ) ) / U11
      U22 = TMP( LOCU22( IPIV ) ) - U12*L21
      XSWAP = XSWPIV( IPIV )
      BSWAP = BSWPIV( IPIV )
      IF( ABS( U22 )<=SMIN ) THEN
         INFO = 1
         U22 = SMIN
      END IF
      IF( BSWAP ) THEN
         TEMP = BTMP( 2 )
         BTMP( 2 ) = BTMP( 1 ) - L21*TEMP
         BTMP( 1 ) = TEMP
      ELSE
         BTMP( 2 ) = BTMP( 2 ) - L21*BTMP( 1 )
      END IF
      SCALE = ONE
      IF( ( TWO*SMLNUM )*ABS( BTMP( 2 ) )>ABS( U22 ) .OR. &
          ( TWO*SMLNUM )*ABS( BTMP( 1 ) )>ABS( U11 ) ) THEN
         SCALE = HALF / MAX( ABS( BTMP( 1 ) ), ABS( BTMP( 2 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
      END IF
      X2( 2 ) = BTMP( 2 ) / U22
      X2( 1 ) = BTMP( 1 ) / U11 - ( U12 / U11 )*X2( 2 )
      IF( XSWAP ) THEN
         TEMP = X2( 2 )
         X2( 2 ) = X2( 1 )
         X2( 1 ) = TEMP
      END IF
      X( 1, 1 ) = X2( 1 )
      IF( N1==1 ) THEN
         X( 1, 2 ) = X2( 2 )
         XNORM = ABS( X( 1, 1 ) ) + ABS( X( 1, 2 ) )
      ELSE
         X( 2, 1 ) = X2( 2 )
         XNORM = MAX( ABS( X( 1, 1 ) ), ABS( X( 2, 1 ) ) )
      END IF
      RETURN
!
!     2 by 2:
!     op[TL11 TL12]*[X11 X12] +ISGN* [X11 X12]*op[TR11 TR12] = [B11 B12]
!       [TL21 TL22] [X21 X22]        [X21 X22]   [TR21 TR22]   [B21 B22]
!
!     Solve equivalent 4 by 4 system using complete pivoting.
!     Set pivots less than SMIN to SMIN.
!
   50 CONTINUE
      SMIN = MAX( ABS( TR( 1, 1 ) ), ABS( TR( 1, 2 ) ), &
             ABS( TR( 2, 1 ) ), ABS( TR( 2, 2 ) ) )
      SMIN = MAX( SMIN, ABS( TL( 1, 1 ) ), ABS( TL( 1, 2 ) ), &
             ABS( TL( 2, 1 ) ), ABS( TL( 2, 2 ) ) )
      SMIN = MAX( EPS*SMIN, SMLNUM )
      BTMP( 1 ) = ZERO
      CALL SCOPY( 16, BTMP, 0, T16, 1 )
      T16( 1, 1 ) = TL( 1, 1 ) + SGN*TR( 1, 1 )
      T16( 2, 2 ) = TL( 2, 2 ) + SGN*TR( 1, 1 )
      T16( 3, 3 ) = TL( 1, 1 ) + SGN*TR( 2, 2 )
      T16( 4, 4 ) = TL( 2, 2 ) + SGN*TR( 2, 2 )
      IF( LTRANL ) THEN
         T16( 1, 2 ) = TL( 2, 1 )
         T16( 2, 1 ) = TL( 1, 2 )
         T16( 3, 4 ) = TL( 2, 1 )
         T16( 4, 3 ) = TL( 1, 2 )
      ELSE
         T16( 1, 2 ) = TL( 1, 2 )
         T16( 2, 1 ) = TL( 2, 1 )
         T16( 3, 4 ) = TL( 1, 2 )
         T16( 4, 3 ) = TL( 2, 1 )
      END IF
      IF( LTRANR ) THEN
         T16( 1, 3 ) = SGN*TR( 1, 2 )
         T16( 2, 4 ) = SGN*TR( 1, 2 )
         T16( 3, 1 ) = SGN*TR( 2, 1 )
         T16( 4, 2 ) = SGN*TR( 2, 1 )
      ELSE
         T16( 1, 3 ) = SGN*TR( 2, 1 )
         T16( 2, 4 ) = SGN*TR( 2, 1 )
         T16( 3, 1 ) = SGN*TR( 1, 2 )
         T16( 4, 2 ) = SGN*TR( 1, 2 )
      END IF
      BTMP( 1 ) = B( 1, 1 )
      BTMP( 2 ) = B( 2, 1 )
      BTMP( 3 ) = B( 1, 2 )
      BTMP( 4 ) = B( 2, 2 )
!
!     Perform elimination
!
      DO 100 I = 1, 3
         XMAX = ZERO
         DO 70 IP = I, 4
            DO 60 JP = I, 4
               IF( ABS( T16( IP, JP ) )>=XMAX ) THEN
                  XMAX = ABS( T16( IP, JP ) )
                  IPSV = IP
                  JPSV = JP
               END IF
   60       CONTINUE
   70    CONTINUE
         IF( IPSV/=I ) THEN
            CALL SSWAP( 4, T16( IPSV, 1 ), 4, T16( I, 1 ), 4 )
            TEMP = BTMP( I )
            BTMP( I ) = BTMP( IPSV )
            BTMP( IPSV ) = TEMP
         END IF
         IF( JPSV/=I ) &
            CALL SSWAP( 4, T16( 1, JPSV ), 1, T16( 1, I ), 1 )
         JPIV( I ) = JPSV
         IF( ABS( T16( I, I ) )<SMIN ) THEN
            INFO = 1
            T16( I, I ) = SMIN
         END IF
         DO 90 J = I + 1, 4
            T16( J, I ) = T16( J, I ) / T16( I, I )
            BTMP( J ) = BTMP( J ) - T16( J, I )*BTMP( I )
            DO 80 K = I + 1, 4
               T16( J, K ) = T16( J, K ) - T16( J, I )*T16( I, K )
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      IF( ABS( T16( 4, 4 ) )<SMIN ) &
         T16( 4, 4 ) = SMIN
      SCALE = ONE
      IF( ( EIGHT*SMLNUM )*ABS( BTMP( 1 ) )>ABS( T16( 1, 1 ) ) .OR. &
          ( EIGHT*SMLNUM )*ABS( BTMP( 2 ) )>ABS( T16( 2, 2 ) ) .OR. &
          ( EIGHT*SMLNUM )*ABS( BTMP( 3 ) )>ABS( T16( 3, 3 ) ) .OR. &
          ( EIGHT*SMLNUM )*ABS( BTMP( 4 ) )>ABS( T16( 4, 4 ) ) ) THEN
         SCALE = ( ONE / EIGHT ) / MAX( ABS( BTMP( 1 ) ), &
                 ABS( BTMP( 2 ) ), ABS( BTMP( 3 ) ), ABS( BTMP( 4 ) ) )
         BTMP( 1 ) = BTMP( 1 )*SCALE
         BTMP( 2 ) = BTMP( 2 )*SCALE
         BTMP( 3 ) = BTMP( 3 )*SCALE
         BTMP( 4 ) = BTMP( 4 )*SCALE
      END IF
      DO 120 I = 1, 4
         K = 5 - I
         TEMP = ONE / T16( K, K )
         TMP( K ) = BTMP( K )*TEMP
         DO 110 J = K + 1, 4
            TMP( K ) = TMP( K ) - ( TEMP*T16( K, J ) )*TMP( J )
  110    CONTINUE
  120 CONTINUE
      DO 130 I = 1, 3
         IF( JPIV( 4-I )/=4-I ) THEN
            TEMP = TMP( 4-I )
            TMP( 4-I ) = TMP( JPIV( 4-I ) )
            TMP( JPIV( 4-I ) ) = TEMP
         END IF
  130 CONTINUE
      X( 1, 1 ) = TMP( 1 )
      X( 2, 1 ) = TMP( 2 )
      X( 1, 2 ) = TMP( 3 )
      X( 2, 2 ) = TMP( 4 )
      XNORM = MAX( ABS( TMP( 1 ) )+ABS( TMP( 3 ) ), &
              ABS( TMP( 2 ) )+ABS( TMP( 4 ) ) )
      RETURN
!
!     End of SLASY2
!
      END
      SUBROUTINE SORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, &
                         WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  SORM2R overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by SGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) REAL array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          SGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) REAL array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by SGEQRF.
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) REAL array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      REAL               AII
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M<0 ) THEN
         INFO = -3
      ELSE IF( N<0 ) THEN
         INFO = -4
      ELSE IF( K<0 .OR. K>NQ ) THEN
         INFO = -5
      ELSE IF( LDA<MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC<MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SORM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M==0 .OR. N==0 .OR. K==0 ) &
         RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) &
           THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i)
!
         AII = A( I, I )
         A( I, I ) = ONE
         CALL SLARF( SIDE, MI, NI, A( I, I ), 1, TAU( I ), C( IC, JC ), &
                     LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of SORM2R
!
      END
      SUBROUTINE SPTTRF( N, D, E, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      REAL               D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  SPTTRF computes the factorization of a real symmetric positive
!  definite tridiagonal matrix A.
!
!  If the subdiagonal elements of A are supplied in the array E, the
!  factorization has the form A = L*D*L**T, where D is diagonal and L
!  is unit lower bidiagonal; if the superdiagonal elements of A are
!  supplied, it has the form A = U**T*D*U, where U is unit upper
!  bidiagonal.  (The two forms are equivalent if A is real.)
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  D       (input/output) REAL array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix
!          A.  On exit, the n diagonal elements of the diagonal matrix
!          D from the L*D*L**T factorization of A.
!
!  E       (input/output) REAL array, dimension (N-1)
!          On entry, the (n-1) off-diagonal elements of the tridiagonal
!          matrix A.
!          On exit, the (n-1) off-diagonal elements of the unit
!          bidiagonal factor L or U from the factorization of A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite; if i < N, the factorization could
!                not be completed, while if i = N, the factorization was
!                completed, but D(N) = 0.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      INTEGER            I
      REAL               DI, EI
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF( N<0 ) THEN
         INFO = -1
         CALL XERBLA( 'SPTTRF', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N==0 ) &
         RETURN
!
!     Compute the L*D*L' (or U'*D*U) factorization of A.
!
      DO 10 I = 1, N - 1
!
!        Drop out of the loop if d(i) <= 0: the matrix is not positive
!        definite.
!
         DI = D( I )
         IF( DI<=ZERO ) &
            GO TO 20
!
!        Solve for e(i) and d(i+1).
!
         EI = E( I )
         E( I ) = EI / DI
         D( I+1 ) = D( I+1 ) - E( I )*EI
   10 CONTINUE
!
!     Check d(n) for positive definiteness.
!
      I = N
      IF( D( I )>ZERO ) &
         GO TO 30
!
   20 CONTINUE
      INFO = I
!
   30 CONTINUE
      RETURN
!
!     End of SPTTRF
!
      END
      SUBROUTINE SPTTRS( N, NRHS, D, E, B, LDB, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL               B( LDB, * ), D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  SPTTRS solves a system of linear equations A * X = B with a
!  symmetric positive definite tridiagonal matrix A using the
!  factorization A = L*D*L**T or A = U**T*D*U computed by SPTTRF.
!  (The two forms are equivalent if A is real.)
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the tridiagonal matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  D       (input) REAL array, dimension (N)
!          The n diagonal elements of the diagonal matrix D from the
!          factorization computed by SPTTRF.
!
!  E       (input) REAL array, dimension (N-1)
!          The (n-1) off-diagonal elements of the unit bidiagonal factor
!          U or L from the factorization computed by SPTTRF.
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            I, J
!     ..
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      INFO = 0
      IF( N<0 ) THEN
         INFO = -1
      ELSE IF( NRHS<0 ) THEN
         INFO = -2
      ELSE IF( LDB<MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SPTTRS', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N==0 ) &
         RETURN
!
!     Solve A * X = B using the factorization A = L*D*L',
!     overwriting each right hand side vector with its solution.
!
      DO 30 J = 1, NRHS
!
!        Solve L * x = b.
!
         DO 10 I = 2, N
            B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   10    CONTINUE
!
!        Solve D * L' * x = b.
!
         B( N, J ) = B( N, J ) / D( N )
         DO 20 I = N - 1, 1, -1
            B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   20    CONTINUE
   30 CONTINUE
!
      RETURN
!
!     End of SPTTRS
!
      END
      SUBROUTINE SSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
!     ..
!     .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  SSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band symmetric matrix can also be found
!  if SSYTRD or SSPTRD or SSBTRD has been used to reduce this matrix to
!  tridiagonal form.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  symmetric matrix.  On entry, Z must contain the
!                  orthogonal matrix used to reduce the original matrix
!                  to tridiagonal form.
!          = 'I':  Compute eigenvalues and eigenvectors of the
!                  tridiagonal matrix.  Z is initialized to the identity
!                  matrix.
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) REAL array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) REAL array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) REAL array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).
!
!  WORK    (workspace) REAL array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero; on exit, D
!                and E contain the elements of a symmetric tridiagonal
!                matrix which is orthogonally similar to the original
!                matrix.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0, &
                         THREE = 3.0E0 )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, &
                         LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1, &
                         NM1, NMAXIT
      REAL               ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2, &
                         S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLAMCH, SLANST, SLAPY2
      EXTERNAL           LSAME, SLAMCH, SLANST, SLAPY2
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAE2, SLAEV2, SLARTG, SLASCL, SLASET, SLASR, &
                         SLASRT, SSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ<0 ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ<1 ) .OR. ( ICOMPZ>0 .AND. LDZ<MAX( 1, &
               N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'SSTEQR', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N==0 ) &
         RETURN
!
      IF( N==1 ) THEN
         IF( ICOMPZ==2 ) &
            Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     Determine the unit roundoff and over/underflow thresholds.
!
      EPS = SLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = SLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.
!
      IF( ICOMPZ==2 ) &
         CALL SLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF( L1>N ) &
         GO TO 160
      IF( L1>1 ) &
         E( L1-1 ) = ZERO
      IF( L1<=NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST==ZERO ) &
               GO TO 30
            IF( TST<=( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+ &
                1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND==L ) &
         GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = SLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM==ZERO ) &
         GO TO 10
      IF( ANORM>SSFMAX ) THEN
         ISCALE = 1
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N, &
                      INFO )
      ELSE IF( ANORM<SSFMIN ) THEN
         ISCALE = 2
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N, &
                      INFO )
         CALL SLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N, &
                      INFO )
      END IF
!
!     Choose between QL and QR iteration
!
      IF( ABS( D( LEND ) )<ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND>L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   40    CONTINUE
         IF( L/=LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST<=( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+ &
                   SAFMIN )GO TO 60
   50       CONTINUE
         END IF
!
         M = LEND
!
   60    CONTINUE
         IF( M<LEND ) &
            E( M ) = ZERO
         P = D( L )
         IF( M==L ) &
            GO TO 80
!
!        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M==L+1 ) THEN
            IF( ICOMPZ>0 ) THEN
               CALL SLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL SLASR( 'R', 'V', 'B', N, 2, WORK( L ), &
                           WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL SLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L<=LEND ) &
               GO TO 40
            GO TO 140
         END IF
!
         IF( JTOT==NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I/=M-1 ) &
               E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ>0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
!
   70    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ>0 ) THEN
            MM = M - L + 1
            CALL SLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ), &
                        Z( 1, L ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
!
!        Eigenvalue found.
!
   80    CONTINUE
         D( L ) = P
!
         L = L + 1
         IF( L<=LEND ) &
            GO TO 40
         GO TO 140
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
   90    CONTINUE
         IF( L/=LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST<=( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+ &
                   SAFMIN )GO TO 110
  100       CONTINUE
         END IF
!
         M = LEND
!
  110    CONTINUE
         IF( M>LEND ) &
            E( M-1 ) = ZERO
         P = D( L )
         IF( M==L ) &
            GO TO 130
!
!        If remaining matrix is 2-by-2, use SLAE2 or SLAEV2
!        to compute its eigensystem.
!
         IF( M==L-1 ) THEN
            IF( ICOMPZ>0 ) THEN
               CALL SLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL SLASR( 'R', 'V', 'F', N, 2, WORK( M ), &
                           WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL SLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L>=LEND ) &
               GO TO 90
            GO TO 140
         END IF
!
         IF( JTOT==NMAXIT ) &
            GO TO 140
         JTOT = JTOT + 1
!
!        Form shift.
!
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = SLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
!
         S = ONE
         C = ONE
         P = ZERO
!
!        Inner loop
!
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL SLARTG( G, F, C, S, R )
            IF( I/=M ) &
               E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
!
!           If eigenvectors are desired, then save rotations.
!
            IF( ICOMPZ>0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
!
  120    CONTINUE
!
!        If eigenvectors are desired, then apply saved rotations.
!
         IF( ICOMPZ>0 ) THEN
            MM = L - M + 1
            CALL SLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ), &
                        Z( 1, M ), LDZ )
         END IF
!
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
!
!        Eigenvalue found.
!
  130    CONTINUE
         D( L ) = P
!
         L = L - 1
         IF( L>=LEND ) &
            GO TO 90
         GO TO 140
!
      END IF
!
!     Undo scaling if necessary
!
  140 CONTINUE
      IF( ISCALE==1 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      ELSE IF( ISCALE==2 ) THEN
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1, &
                      D( LSV ), N, INFO )
         CALL SLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ), &
                      N, INFO )
      END IF
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT<NMAXIT ) &
         GO TO 10
      DO 150 I = 1, N - 1
         IF( E( I )/=ZERO ) &
            INFO = INFO + 1
  150 CONTINUE
      GO TO 190
!
!     Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      IF( ICOMPZ==0 ) THEN
!
!        Use Quick Sort
!
         CALL SLASRT( 'I', N, D, INFO )
!
      ELSE
!
!        Use Selection Sort to minimize swaps of eigenvectors
!
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J )<P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K/=I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL SSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
!
  190 CONTINUE
      RETURN
!
!     End of SSTEQR
!
      END
      SUBROUTINE STREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR, &
                         LDVR, MM, M, WORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      REAL               T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  STREVC computes some or all of the right and/or left eigenvectors of
!  a real upper quasi-triangular matrix T.
!
!  The right eigenvector x and the left eigenvector y of T corresponding
!  to an eigenvalue w are defined by:
!
!               T*x = w*x,     y'*T = w*y'
!
!  where y' denotes the conjugate transpose of the vector y.
!
!  If all eigenvectors are requested, the routine may either return the
!  matrices X and/or Y of right or left eigenvectors of T, or the
!  products Q*X and/or Q*Y, where Q is an input orthogonal
!  matrix. If T was obtained from the real-Schur factorization of an
!  original matrix A = Q*T*Q', then Q*X and Q*Y are the matrices of
!  right or left eigenvectors of A.
!
!  T must be in Schur canonical form (as returned by SHSEQR), that is,
!  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!  2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.  Corresponding to each 2-by-2
!  diagonal block is a complex conjugate pair of eigenvalues and
!  eigenvectors; only one eigenvector of the pair is computed, namely
!  the one corresponding to the eigenvalue with positive imaginary part.

!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'R':  compute right eigenvectors only;
!          = 'L':  compute left eigenvectors only;
!          = 'B':  compute both right and left eigenvectors.
!
!  HOWMNY  (input) CHARACTER*1
!          = 'A':  compute all right and/or left eigenvectors;
!          = 'B':  compute all right and/or left eigenvectors,
!                  and backtransform them using the input matrices
!                  supplied in VR and/or VL;
!          = 'S':  compute selected right and/or left eigenvectors,
!                  specified by the logical array SELECT.
!
!  SELECT  (input/output) LOGICAL array, dimension (N)
!          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!          computed.
!          If HOWMNY = 'A' or 'B', SELECT is not referenced.
!          To select the real eigenvector corresponding to a real
!          eigenvalue w(j), SELECT(j) must be set to .TRUE..  To select
!          the complex eigenvector corresponding to a complex conjugate
!          pair w(j) and w(j+1), either SELECT(j) or SELECT(j+1) must be
!          set to .TRUE.; then on exit SELECT(j) is .TRUE. and
!          SELECT(j+1) is .FALSE..
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input) REAL array, dimension (LDT,N)
!          The upper quasi-triangular matrix T in Schur canonical form.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  VL      (input/output) REAL array, dimension (LDVL,MM)
!          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!          of Schur vectors returned by SHSEQR).
!          On exit, if SIDE = 'L' or 'B', VL contains:
!          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*Y;
!          if HOWMNY = 'S', the left eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VL, in the same order as their
!                           eigenvalues.
!          A complex eigenvector corresponding to a complex eigenvalue
!          is stored in two consecutive columns, the first holding the
!          real part, and the second the imaginary part.
!          If SIDE = 'R', VL is not referenced.
!
!  LDVL    (input) INTEGER
!          The leading dimension of the array VL.  LDVL >= max(1,N) if
!          SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
!
!  VR      (input/output) REAL array, dimension (LDVR,MM)
!          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!          of Schur vectors returned by SHSEQR).
!          On exit, if SIDE = 'R' or 'B', VR contains:
!          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!          if HOWMNY = 'B', the matrix Q*X;
!          if HOWMNY = 'S', the right eigenvectors of T specified by
!                           SELECT, stored consecutively in the columns
!                           of VR, in the same order as their
!                           eigenvalues.
!          A complex eigenvector corresponding to a complex eigenvalue
!          is stored in two consecutive columns, the first holding the
!          real part and the second the imaginary part.
!          If SIDE = 'L', VR is not referenced.
!
!  LDVR    (input) INTEGER
!          The leading dimension of the array VR.  LDVR >= max(1,N) if
!          SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
!
!  MM      (input) INTEGER
!          The number of columns in the arrays VL and/or VR. MM >= M.
!
!  M       (output) INTEGER
!          The number of columns in the arrays VL and/or VR actually
!          used to store the eigenvectors.
!          If HOWMNY = 'A' or 'B', M is set to N.
!          Each selected real eigenvector occupies one column and each
!          selected complex eigenvector occupies two columns.
!
!  WORK    (workspace) REAL array, dimension (3*N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  The algorithm used in this program is basically backward (forward)
!  substitution, with scaling to make the the code robust against
!  possible overflow.
!
!  Each eigenvector is normalized so that the element of largest
!  magnitude has magnitude 1; here the magnitude of a complex number
!  (x,y) is taken to be |x| + |y|.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            ALLV, BOTHV, LEFTV, OVER, PAIR, RIGHTV, SOMEV
      INTEGER            I, IERR, II, IP, IS, J, J1, J2, JNXT, K, KI, N2
      REAL               BETA, BIGNUM, EMAX, OVFL, REC, REMAX, SCALE, &
                         SMIN, SMLNUM, ULP, UNFL, VCRIT, VMAX, WI, WR, &
                         XNORM
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ISAMAX
      REAL               SDOT, SLAMCH
      EXTERNAL           LSAME, ISAMAX, SDOT, SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SGEMV, SLABAD, SLALN2, SSCAL, &
                         XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Local Arrays ..
      REAL               X( 2, 2 )
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      BOTHV = LSAME( SIDE, 'B' )
      RIGHTV = LSAME( SIDE, 'R' ) .OR. BOTHV
      LEFTV = LSAME( SIDE, 'L' ) .OR. BOTHV
!
      ALLV = LSAME( HOWMNY, 'A' )
      OVER = LSAME( HOWMNY, 'B' ) .OR. LSAME( HOWMNY, 'O' )
      SOMEV = LSAME( HOWMNY, 'S' )
!
      INFO = 0
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.ALLV .AND. .NOT.OVER .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N<0 ) THEN
         INFO = -4
      ELSE IF( LDT<MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL<1 .OR. ( LEFTV .AND. LDVL<N ) ) THEN
         INFO = -8
      ELSE IF( LDVR<1 .OR. ( RIGHTV .AND. LDVR<N ) ) THEN
         INFO = -10
      ELSE
!
!        Set M to the number of columns required to store the selected
!        eigenvectors, standardize the array SELECT if necessary, and
!        test MM.
!
         IF( SOMEV ) THEN
            M = 0
            PAIR = .FALSE.
            DO 10 J = 1, N
               IF( PAIR ) THEN
                  PAIR = .FALSE.
                  SELECT( J ) = .FALSE.
               ELSE
                  IF( J<N ) THEN
                     IF( T( J+1, J )==ZERO ) THEN
                        IF( SELECT( J ) ) &
                           M = M + 1
                     ELSE
                        PAIR = .TRUE.
                        IF( SELECT( J ) .OR. SELECT( J+1 ) ) THEN
                           SELECT( J ) = .TRUE.
                           M = M + 2
                        END IF
                     END IF
                  ELSE
                     IF( SELECT( N ) ) &
                        M = M + 1
                  END IF
               END IF
   10       CONTINUE
         ELSE
            M = N
         END IF
!
         IF( MM<M ) THEN
            INFO = -11
         END IF
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'STREVC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( N==0 ) &
         RETURN
!
!     Set the constants to control overflow.
!
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )
      BIGNUM = ( ONE-ULP ) / SMLNUM
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      WORK( 1 ) = ZERO
      DO 30 J = 2, N
         WORK( J ) = ZERO
         DO 20 I = 1, J - 1
            WORK( J ) = WORK( J ) + ABS( T( I, J ) )
   20    CONTINUE
   30 CONTINUE
!
!     Index IP is used to specify the real or complex eigenvalue:
!       IP = 0, real eigenvalue,
!            1, first of conjugate complex pair: (wr,wi)
!           -1, second of conjugate complex pair: (wr,wi)
!
      N2 = 2*N
!
      IF( RIGHTV ) THEN
!
!        Compute right eigenvectors.
!
         IP = 0
         IS = M
         DO 140 KI = N, 1, -1
!
            IF( IP==1 ) &
               GO TO 130
            IF( KI==1 ) &
               GO TO 40
            IF( T( KI, KI-1 )==ZERO ) &
               GO TO 40
            IP = -1
!
   40       CONTINUE
            IF( SOMEV ) THEN
               IF( IP==0 ) THEN
                  IF( .NOT.SELECT( KI ) ) &
                     GO TO 130
               ELSE
                  IF( .NOT.SELECT( KI-1 ) ) &
                     GO TO 130
               END IF
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            WR = T( KI, KI )
            WI = ZERO
            IF( IP/=0 ) &
               WI = SQRT( ABS( T( KI, KI-1 ) ) )* &
                    SQRT( ABS( T( KI-1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
!
            IF( IP==0 ) THEN
!
!              Real right eigenvector
!
               WORK( KI+N ) = ONE
!
!              Form right-hand side
!
               DO 50 K = 1, KI - 1
                  WORK( K+N ) = -T( K, KI )
   50          CONTINUE
!
!              Solve the upper quasi-triangular system:
!                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
!
               JNXT = KI - 1
               DO 60 J = KI - 1, 1, -1
                  IF( J>JNXT ) &
                     GO TO 60
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J>1 ) THEN
                     IF( T( J, J-1 )/=ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
!
                  IF( J1==J2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL SLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale X(1,1) to avoid overflow when updating
!                    the right-hand side.
!
                     IF( XNORM>ONE ) THEN
                        IF( WORK( J )>BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) &
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
!
!                    Update right-hand side
!
                     CALL SAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, &
                                  T( J-1, J-1 ), LDT, ONE, ONE, &
                                  WORK( J-1+N ), N, WR, ZERO, X, 2, &
                                  SCALE, XNORM, IERR )
!
!                    Scale X(1,1) and X(2,1) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( XNORM>ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA>BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 2, 1 ) = X( 2, 1 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) &
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
!
!                    Update right-hand side
!
                     CALL SAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                  END IF
   60          CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( KI, WORK( 1+N ), 1, VR( 1, IS ), 1 )
!
                  II = ISAMAX( KI, VR( 1, IS ), 1 )
                  REMAX = ONE / ABS( VR( II, IS ) )
                  CALL SSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
                  DO 70 K = KI + 1, N
                     VR( K, IS ) = ZERO
   70             CONTINUE
               ELSE
                  IF( KI>1 ) &
                     CALL SGEMV( 'N', N, KI-1, ONE, VR, LDVR, &
                                 WORK( 1+N ), 1, WORK( KI+N ), &
                                 VR( 1, KI ), 1 )
!
                  II = ISAMAX( N, VR( 1, KI ), 1 )
                  REMAX = ONE / ABS( VR( II, KI ) )
                  CALL SSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
!
            ELSE
!
!              Complex right eigenvector.
!
!              Initial solve
!                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
!                [ (T(KI,KI-1)   T(KI,KI)   )               ]
!
               IF( ABS( T( KI-1, KI ) )>=ABS( T( KI, KI-1 ) ) ) THEN
                  WORK( KI-1+N ) = ONE
                  WORK( KI+N2 ) = WI / T( KI-1, KI )
               ELSE
                  WORK( KI-1+N ) = -WI / T( KI, KI-1 )
                  WORK( KI+N2 ) = ONE
               END IF
               WORK( KI+N ) = ZERO
               WORK( KI-1+N2 ) = ZERO
!
!              Form right-hand side
!
               DO 80 K = 1, KI - 2
                  WORK( K+N ) = -WORK( KI-1+N )*T( K, KI-1 )
                  WORK( K+N2 ) = -WORK( KI+N2 )*T( K, KI )
   80          CONTINUE
!
!              Solve upper quasi-triangular system:
!              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
!
               JNXT = KI - 2
               DO 90 J = KI - 2, 1, -1
                  IF( J>JNXT ) &
                     GO TO 90
                  J1 = J
                  J2 = J
                  JNXT = J - 1
                  IF( J>1 ) THEN
                     IF( T( J, J-1 )/=ZERO ) THEN
                        J1 = J - 1
                        JNXT = J - 2
                     END IF
                  END IF
!
                  IF( J1==J2 ) THEN
!
!                    1-by-1 diagonal block
!
                     CALL SLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, WI, &
                                  X, 2, SCALE, XNORM, IERR )
!
!                    Scale X(1,1) and X(1,2) to avoid overflow when
!                    updating the right-hand side.
!
                     IF( XNORM>ONE ) THEN
                        IF( WORK( J )>BIGNUM / XNORM ) THEN
                           X( 1, 1 ) = X( 1, 1 ) / XNORM
                           X( 1, 2 ) = X( 1, 2 ) / XNORM
                           SCALE = SCALE / XNORM
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) THEN
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL SSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
!
!                    Update the right-hand side
!
                     CALL SAXPY( J-1, -X( 1, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-1, -X( 1, 2 ), T( 1, J ), 1, &
                                 WORK( 1+N2 ), 1 )
!
                  ELSE
!
!                    2-by-2 diagonal block
!
                     CALL SLALN2( .FALSE., 2, 2, SMIN, ONE, &
                                  T( J-1, J-1 ), LDT, ONE, ONE, &
                                  WORK( J-1+N ), N, WR, WI, X, 2, SCALE, &
                                  XNORM, IERR )
!
!                    Scale X to avoid overflow when updating
!                    the right-hand side.
!
                     IF( XNORM>ONE ) THEN
                        BETA = MAX( WORK( J-1 ), WORK( J ) )
                        IF( BETA>BIGNUM / XNORM ) THEN
                           REC = ONE / XNORM
                           X( 1, 1 ) = X( 1, 1 )*REC
                           X( 1, 2 ) = X( 1, 2 )*REC
                           X( 2, 1 ) = X( 2, 1 )*REC
                           X( 2, 2 ) = X( 2, 2 )*REC
                           SCALE = SCALE*REC
                        END IF
                     END IF
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) THEN
                        CALL SSCAL( KI, SCALE, WORK( 1+N ), 1 )
                        CALL SSCAL( KI, SCALE, WORK( 1+N2 ), 1 )
                     END IF
                     WORK( J-1+N ) = X( 1, 1 )
                     WORK( J+N ) = X( 2, 1 )
                     WORK( J-1+N2 ) = X( 1, 2 )
                     WORK( J+N2 ) = X( 2, 2 )
!
!                    Update the right-hand side
!
                     CALL SAXPY( J-2, -X( 1, 1 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-2, -X( 2, 1 ), T( 1, J ), 1, &
                                 WORK( 1+N ), 1 )
                     CALL SAXPY( J-2, -X( 1, 2 ), T( 1, J-1 ), 1, &
                                 WORK( 1+N2 ), 1 )
                     CALL SAXPY( J-2, -X( 2, 2 ), T( 1, J ), 1, &
                                 WORK( 1+N2 ), 1 )
                  END IF
   90          CONTINUE
!
!              Copy the vector x or Q*x to VR and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( KI, WORK( 1+N ), 1, VR( 1, IS-1 ), 1 )
                  CALL SCOPY( KI, WORK( 1+N2 ), 1, VR( 1, IS ), 1 )
!
                  EMAX = ZERO
                  DO 100 K = 1, KI
                     EMAX = MAX( EMAX, ABS( VR( K, IS-1 ) )+ &
                            ABS( VR( K, IS ) ) )
  100             CONTINUE
!
                  REMAX = ONE / EMAX
                  CALL SSCAL( KI, REMAX, VR( 1, IS-1 ), 1 )
                  CALL SSCAL( KI, REMAX, VR( 1, IS ), 1 )
!
                  DO 110 K = KI + 1, N
                     VR( K, IS-1 ) = ZERO
                     VR( K, IS ) = ZERO
  110             CONTINUE
!
               ELSE
!
                  IF( KI>2 ) THEN
                     CALL SGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                 WORK( 1+N ), 1, WORK( KI-1+N ), &
                                 VR( 1, KI-1 ), 1 )
                     CALL SGEMV( 'N', N, KI-2, ONE, VR, LDVR, &
                                 WORK( 1+N2 ), 1, WORK( KI+N2 ), &
                                 VR( 1, KI ), 1 )
                  ELSE
                     CALL SSCAL( N, WORK( KI-1+N ), VR( 1, KI-1 ), 1 )
                     CALL SSCAL( N, WORK( KI+N2 ), VR( 1, KI ), 1 )
                  END IF
!
                  EMAX = ZERO
                  DO 120 K = 1, N
                     EMAX = MAX( EMAX, ABS( VR( K, KI-1 ) )+ &
                            ABS( VR( K, KI ) ) )
  120             CONTINUE
                  REMAX = ONE / EMAX
                  CALL SSCAL( N, REMAX, VR( 1, KI-1 ), 1 )
                  CALL SSCAL( N, REMAX, VR( 1, KI ), 1 )
               END IF
            END IF
!
            IS = IS - 1
            IF( IP/=0 ) &
               IS = IS - 1
  130       CONTINUE
            IF( IP==1 ) &
               IP = 0
            IF( IP==-1 ) &
               IP = 1
  140    CONTINUE
      END IF
!
      IF( LEFTV ) THEN
!
!        Compute left eigenvectors.
!
         IP = 0
         IS = 1
         DO 260 KI = 1, N
!
            IF( IP==-1 ) &
               GO TO 250
            IF( KI==N ) &
               GO TO 150
            IF( T( KI+1, KI )==ZERO ) &
               GO TO 150
            IP = 1
!
  150       CONTINUE
            IF( SOMEV ) THEN
               IF( .NOT.SELECT( KI ) ) &
                  GO TO 250
            END IF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
            WR = T( KI, KI )
            WI = ZERO
            IF( IP/=0 ) &
               WI = SQRT( ABS( T( KI, KI+1 ) ) )* &
                    SQRT( ABS( T( KI+1, KI ) ) )
            SMIN = MAX( ULP*( ABS( WR )+ABS( WI ) ), SMLNUM )
!
            IF( IP==0 ) THEN
!
!              Real left eigenvector.
!
               WORK( KI+N ) = ONE
!
!              Form right-hand side
!
               DO 160 K = KI + 1, N
                  WORK( K+N ) = -T( KI, K )
  160          CONTINUE
!
!              Solve the quasi-triangular system:
!                 (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
!
               VMAX = ONE
               VCRIT = BIGNUM
!
               JNXT = KI + 1
               DO 170 J = KI + 1, N
                  IF( J<JNXT ) &
                     GO TO 170
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J<N ) THEN
                     IF( T( J+1, J )/=ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
!
                  IF( J1==J2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     IF( WORK( J )>VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-1, T( KI+1, J ), 1, &
                                   WORK( KI+1+N ), 1 )
!
!                    Solve (T(J,J)-WR)'*X = WORK
!
                     CALL SLALN2( .FALSE., 1, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) &
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     VMAX = MAX( ABS( WORK( J+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA>VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-1, T( KI+1, J ), 1, &
                                   WORK( KI+1+N ), 1 )
!
                     WORK( J+1+N ) = WORK( J+1+N ) - &
                                     SDOT( J-KI-1, T( KI+1, J+1 ), 1, &
                                     WORK( KI+1+N ), 1 )
!
!                    Solve
!                      [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
!                      [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
!
                     CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  ZERO, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) &
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+1+N ) = X( 2, 1 )
!
                     VMAX = MAX( ABS( WORK( J+N ) ), &
                            ABS( WORK( J+1+N ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  END IF
  170          CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
!
                  II = ISAMAX( N-KI+1, VL( KI, IS ), 1 ) + KI - 1
                  REMAX = ONE / ABS( VL( II, IS ) )
                  CALL SSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
!
                  DO 180 K = 1, KI - 1
                     VL( K, IS ) = ZERO
  180             CONTINUE
!
               ELSE
!
                  IF( KI<N ) &
                     CALL SGEMV( 'N', N, N-KI, ONE, VL( 1, KI+1 ), LDVL, &
                                 WORK( KI+1+N ), 1, WORK( KI+N ), &
                                 VL( 1, KI ), 1 )
!
                  II = ISAMAX( N, VL( 1, KI ), 1 )
                  REMAX = ONE / ABS( VL( II, KI ) )
                  CALL SSCAL( N, REMAX, VL( 1, KI ), 1 )
!
               END IF
!
            ELSE
!
!              Complex left eigenvector.
!
!               Initial solve:
!                 ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
!                 ((T(KI+1,KI) T(KI+1,KI+1))                )
!
               IF( ABS( T( KI, KI+1 ) )>=ABS( T( KI+1, KI ) ) ) THEN
                  WORK( KI+N ) = WI / T( KI, KI+1 )
                  WORK( KI+1+N2 ) = ONE
               ELSE
                  WORK( KI+N ) = ONE
                  WORK( KI+1+N2 ) = -WI / T( KI+1, KI )
               END IF
               WORK( KI+1+N ) = ZERO
               WORK( KI+N2 ) = ZERO
!
!              Form right-hand side
!
               DO 190 K = KI + 2, N
                  WORK( K+N ) = -WORK( KI+N )*T( KI, K )
                  WORK( K+N2 ) = -WORK( KI+1+N2 )*T( KI+1, K )
  190          CONTINUE
!
!              Solve complex quasi-triangular system:
!              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
!
               VMAX = ONE
               VCRIT = BIGNUM
!
               JNXT = KI + 2
               DO 200 J = KI + 2, N
                  IF( J<JNXT ) &
                     GO TO 200
                  J1 = J
                  J2 = J
                  JNXT = J + 1
                  IF( J<N ) THEN
                     IF( T( J+1, J )/=ZERO ) THEN
                        J2 = J + 1
                        JNXT = J + 2
                     END IF
                  END IF
!
                  IF( J1==J2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when
!                    forming the right-hand side elements.
!
                     IF( WORK( J )>VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-2, T( KI+2, J ), 1, &
                                   WORK( KI+2+N ), 1 )
                     WORK( J+N2 ) = WORK( J+N2 ) - &
                                    SDOT( J-KI-2, T( KI+2, J ), 1, &
                                    WORK( KI+2+N2 ), 1 )
!
!                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
!
                     CALL SLALN2( .FALSE., 1, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  -WI, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) THEN
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     VMAX = MAX( ABS( WORK( J+N ) ), &
                            ABS( WORK( J+N2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side elements.
!
                     BETA = MAX( WORK( J ), WORK( J+1 ) )
                     IF( BETA>VCRIT ) THEN
                        REC = ONE / VMAX
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, REC, WORK( KI+N2 ), 1 )
                        VMAX = ONE
                        VCRIT = BIGNUM
                     END IF
!
                     WORK( J+N ) = WORK( J+N ) - &
                                   SDOT( J-KI-2, T( KI+2, J ), 1, &
                                   WORK( KI+2+N ), 1 )
!
                     WORK( J+N2 ) = WORK( J+N2 ) - &
                                    SDOT( J-KI-2, T( KI+2, J ), 1, &
                                    WORK( KI+2+N2 ), 1 )
!
                     WORK( J+1+N ) = WORK( J+1+N ) - &
                                     SDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                     WORK( KI+2+N ), 1 )
!
                     WORK( J+1+N2 ) = WORK( J+1+N2 ) - &
                                      SDOT( J-KI-2, T( KI+2, J+1 ), 1, &
                                      WORK( KI+2+N2 ), 1 )
!
!                    Solve 2-by-2 complex linear equation
!                      ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
!                      ([T(j+1,j) T(j+1,j+1)]             )
!
                     CALL SLALN2( .TRUE., 2, 2, SMIN, ONE, T( J, J ), &
                                  LDT, ONE, ONE, WORK( J+N ), N, WR, &
                                  -WI, X, 2, SCALE, XNORM, IERR )
!
!                    Scale if necessary
!
                     IF( SCALE/=ONE ) THEN
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N ), 1 )
                        CALL SSCAL( N-KI+1, SCALE, WORK( KI+N2 ), 1 )
                     END IF
                     WORK( J+N ) = X( 1, 1 )
                     WORK( J+N2 ) = X( 1, 2 )
                     WORK( J+1+N ) = X( 2, 1 )
                     WORK( J+1+N2 ) = X( 2, 2 )
                     VMAX = MAX( ABS( X( 1, 1 ) ), ABS( X( 1, 2 ) ), &
                            ABS( X( 2, 1 ) ), ABS( X( 2, 2 ) ), VMAX )
                     VCRIT = BIGNUM / VMAX
!
                  END IF
  200          CONTINUE
!
!              Copy the vector x or Q*x to VL and normalize.
!
  210          CONTINUE
               IF( .NOT.OVER ) THEN
                  CALL SCOPY( N-KI+1, WORK( KI+N ), 1, VL( KI, IS ), 1 )
                  CALL SCOPY( N-KI+1, WORK( KI+N2 ), 1, VL( KI, IS+1 ), &
                              1 )
!
                  EMAX = ZERO
                  DO 220 K = KI, N
                     EMAX = MAX( EMAX, ABS( VL( K, IS ) )+ &
                            ABS( VL( K, IS+1 ) ) )
  220             CONTINUE
                  REMAX = ONE / EMAX
                  CALL SSCAL( N-KI+1, REMAX, VL( KI, IS ), 1 )
                  CALL SSCAL( N-KI+1, REMAX, VL( KI, IS+1 ), 1 )
!
                  DO 230 K = 1, KI - 1
                     VL( K, IS ) = ZERO
                     VL( K, IS+1 ) = ZERO
  230             CONTINUE
               ELSE
                  IF( KI<N-1 ) THEN
                     CALL SGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ), &
                                 LDVL, WORK( KI+2+N ), 1, WORK( KI+N ), &
                                 VL( 1, KI ), 1 )
                     CALL SGEMV( 'N', N, N-KI-1, ONE, VL( 1, KI+2 ), &
                                 LDVL, WORK( KI+2+N2 ), 1, &
                                 WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  ELSE
                     CALL SSCAL( N, WORK( KI+N ), VL( 1, KI ), 1 )
                     CALL SSCAL( N, WORK( KI+1+N2 ), VL( 1, KI+1 ), 1 )
                  END IF
!
                  EMAX = ZERO
                  DO 240 K = 1, N
                     EMAX = MAX( EMAX, ABS( VL( K, KI ) )+ &
                            ABS( VL( K, KI+1 ) ) )
  240             CONTINUE
                  REMAX = ONE / EMAX
                  CALL SSCAL( N, REMAX, VL( 1, KI ), 1 )
                  CALL SSCAL( N, REMAX, VL( 1, KI+1 ), 1 )
!
               END IF
!
            END IF
!
            IS = IS + 1
            IF( IP/=0 ) &
               IS = IS + 1
  250       CONTINUE
            IF( IP==-1 ) &
               IP = 0
            IF( IP==1 ) &
               IP = -1
!
  260    CONTINUE
!
      END IF
!
      RETURN
!
!     End of STREVC
!
      END
      SUBROUTINE STREXC( COMPQ, N, T, LDT, Q, LDQ, IFST, ILST, WORK, &
                         INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ
      INTEGER            IFST, ILST, INFO, LDQ, LDT, N
!     ..
!     .. Array Arguments ..
      REAL               Q( LDQ, * ), T( LDT, * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  STREXC reorders the real Schur factorization of a real matrix
!  A = Q*T*Q**T, so that the diagonal block of T with row index IFST is
!  moved to row ILST.
!
!  The real Schur form T is reordered by an orthogonal similarity
!  transformation Z**T*T*Z, and optionally the matrix Q of Schur vectors
!  is updated by postmultiplying it with Z.
!
!  T must be in Schur canonical form (as returned by SHSEQR), that is,
!  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!  2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  COMPQ   (input) CHARACTER*1
!          = 'V':  update the matrix Q of Schur vectors;
!          = 'N':  do not update Q.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) REAL array, dimension (LDT,N)
!          On entry, the upper quasi-triangular matrix T, in Schur
!          Schur canonical form.
!          On exit, the reordered upper quasi-triangular matrix, again
!          in Schur canonical form.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  Q       (input/output) REAL array, dimension (LDQ,N)
!          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!          orthogonal transformation matrix Z which reorders T.
!          If COMPQ = 'N', Q is not referenced.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q.  LDQ >= max(1,N).
!
!  IFST    (input/output) INTEGER
!  ILST    (input/output) INTEGER
!          Specify the reordering of the diagonal blocks of T.
!          The block with row index IFST is moved to row ILST, by a
!          sequence of transpositions between adjacent blocks.
!          On exit, if IFST pointed on entry to the second row of a
!          2-by-2 block, it is changed to point to the first row; ILST
!          always points to the first row of the block in its final
!          position (which may differ from its input value by +1 or -1).
!          1 <= IFST <= N; 1 <= ILST <= N.
!
!  WORK    (workspace) REAL array, dimension (N)
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          = 1:  two adjacent blocks were too close to swap (the problem
!                is very ill-conditioned); T may have been partially
!                reordered, and ILST points to the first row of the
!                current position of the block being moved.
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            WANTQ
      INTEGER            HERE, NBF, NBL, NBNEXT
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLAEXC, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input arguments.
!
      INFO = 0
      WANTQ = LSAME( COMPQ, 'V' )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( COMPQ, 'N' ) ) THEN
         INFO = -1
      ELSE IF( N<0 ) THEN
         INFO = -2
      ELSE IF( LDT<MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDQ<1 .OR. ( WANTQ .AND. LDQ<MAX( 1, N ) ) ) THEN
         INFO = -6
      ELSE IF( IFST<1 .OR. IFST>N ) THEN
         INFO = -7
      ELSE IF( ILST<1 .OR. ILST>N ) THEN
         INFO = -8
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'STREXC', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N<=1 ) &
         RETURN
!
!     Determine the first row of specified block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF( IFST>1 ) THEN
         IF( T( IFST, IFST-1 )/=ZERO ) &
            IFST = IFST - 1
      END IF
      NBF = 1
      IF( IFST<N ) THEN
         IF( T( IFST+1, IFST )/=ZERO ) &
            NBF = 2
      END IF
!
!     Determine the first row of the final block
!     and find out it is 1 by 1 or 2 by 2.
!
      IF( ILST>1 ) THEN
         IF( T( ILST, ILST-1 )/=ZERO ) &
            ILST = ILST - 1
      END IF
      NBL = 1
      IF( ILST<N ) THEN
         IF( T( ILST+1, ILST )/=ZERO ) &
            NBL = 2
      END IF
!
      IF( IFST==ILST ) &
         RETURN
!
      IF( IFST<ILST ) THEN
!
!        Update ILST
!
         IF( NBF==2 .AND. NBL==1 ) &
            ILST = ILST - 1
         IF( NBF==1 .AND. NBL==2 ) &
            ILST = ILST + 1
!
         HERE = IFST
!
   10    CONTINUE
!
!        Swap block with next one below
!
         IF( NBF==1 .OR. NBF==2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
            NBNEXT = 1
            IF( HERE+NBF+1<=N ) THEN
               IF( T( HERE+NBF+1, HERE+NBF )/=ZERO ) &
                  NBNEXT = 2
            END IF
            CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBF, NBNEXT, &
                         WORK, INFO )
            IF( INFO/=0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE + NBNEXT
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            IF( NBF==2 ) THEN
               IF( T( HERE+1, HERE )==ZERO ) &
                  NBF = 3
            END IF
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            NBNEXT = 1
            IF( HERE+3<=N ) THEN
               IF( T( HERE+3, HERE+2 )/=ZERO ) &
                  NBNEXT = 2
            END IF
            CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, NBNEXT, &
                         WORK, INFO )
            IF( INFO/=0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT==1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, NBNEXT, &
                            WORK, INFO )
               HERE = HERE + 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               IF( T( HERE+2, HERE+1 )==ZERO ) &
                  NBNEXT = 1
               IF( NBNEXT==2 ) THEN
!
!                 2 by 2 Block did not split
!
                  CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, &
                               NBNEXT, WORK, INFO )
                  IF( INFO/=0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE + 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                               WORK, INFO )
                  CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE+1, 1, 1, &
                               WORK, INFO )
                  HERE = HERE + 2
               END IF
            END IF
         END IF
         IF( HERE<ILST ) &
            GO TO 10
!
      ELSE
!
         HERE = IFST
   20    CONTINUE
!
!        Swap block with next one above
!
         IF( NBF==1 .OR. NBF==2 ) THEN
!
!           Current block either 1 by 1 or 2 by 2
!
            NBNEXT = 1
            IF( HERE>=3 ) THEN
               IF( T( HERE-1, HERE-2 )/=ZERO ) &
                  NBNEXT = 2
            END IF
            CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                         NBF, WORK, INFO )
            IF( INFO/=0 ) THEN
               ILST = HERE
               RETURN
            END IF
            HERE = HERE - NBNEXT
!
!           Test if 2 by 2 block breaks into two 1 by 1 blocks
!
            IF( NBF==2 ) THEN
               IF( T( HERE+1, HERE )==ZERO ) &
                  NBF = 3
            END IF
!
         ELSE
!
!           Current block consists of two 1 by 1 blocks each of which
!           must be swapped individually
!
            NBNEXT = 1
            IF( HERE>=3 ) THEN
               IF( T( HERE-1, HERE-2 )/=ZERO ) &
                  NBNEXT = 2
            END IF
            CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-NBNEXT, NBNEXT, &
                         1, WORK, INFO )
            IF( INFO/=0 ) THEN
               ILST = HERE
               RETURN
            END IF
            IF( NBNEXT==1 ) THEN
!
!              Swap two 1 by 1 blocks, no problems possible
!
               CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, NBNEXT, 1, &
                            WORK, INFO )
               HERE = HERE - 1
            ELSE
!
!              Recompute NBNEXT in case 2 by 2 split
!
               IF( T( HERE, HERE-1 )==ZERO ) &
                  NBNEXT = 1
               IF( NBNEXT==2 ) THEN
!
!                 2 by 2 Block did not split
!
                  CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 2, 1, &
                               WORK, INFO )
                  IF( INFO/=0 ) THEN
                     ILST = HERE
                     RETURN
                  END IF
                  HERE = HERE - 2
               ELSE
!
!                 2 by 2 Block did split
!
                  CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE, 1, 1, &
                               WORK, INFO )
                  CALL SLAEXC( WANTQ, N, T, LDT, Q, LDQ, HERE-1, 1, 1, &
                               WORK, INFO )
                  HERE = HERE - 2
               END IF
            END IF
         END IF
         IF( HERE>ILST ) &
            GO TO 20
      END IF
      ILST = HERE
!
      RETURN
!
!     End of STREXC
!
      END
      SUBROUTINE STRSEN( JOB, COMPQ, SELECT, N, T, LDT, Q, LDQ, WR, WI, &
                         M, S, SEP, WORK, LWORK, IWORK, LIWORK, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER          COMPQ, JOB
      INTEGER            INFO, LDQ, LDT, LIWORK, LWORK, M, N
      REAL               S, SEP
!     ..
!     .. Array Arguments ..
      LOGICAL            SELECT( * )
      INTEGER            IWORK( * )
      REAL               Q( LDQ, * ), T( LDT, * ), WI( * ), WORK( * ), &
                         WR( * )
!     ..
!
!  Purpose
!  =======
!
!  STRSEN reorders the real Schur factorization of a real matrix
!  A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
!  the leading diagonal blocks of the upper quasi-triangular matrix T,
!  and the leading columns of Q form an orthonormal basis of the
!  corresponding right invariant subspace.
!
!  Optionally the routine computes the reciprocal condition numbers of
!  the cluster of eigenvalues and/or the invariant subspace.
!
!  T must be in Schur canonical form (as returned by SHSEQR), that is,
!  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!  2-by-2 diagonal block has its diagonal elemnts equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  JOB     (input) CHARACTER*1
!          Specifies whether condition numbers are required for the
!          cluster of eigenvalues (S) or the invariant subspace (SEP):
!          = 'N': none;
!          = 'E': for eigenvalues only (S);
!          = 'V': for invariant subspace only (SEP);
!          = 'B': for both eigenvalues and invariant subspace (S and
!                 SEP).
!
!  COMPQ   (input) CHARACTER*1
!          = 'V': update the matrix Q of Schur vectors;
!          = 'N': do not update Q.
!
!  SELECT  (input) LOGICAL array, dimension (N)
!          SELECT specifies the eigenvalues in the selected cluster. To
!          select a real eigenvalue w(j), SELECT(j) must be set to
!          .TRUE.. To select a complex conjugate pair of eigenvalues
!          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
!          either SELECT(j) or SELECT(j+1) or both must be set to
!          .TRUE.; a complex conjugate pair of eigenvalues must be
!          either both included in the cluster or both excluded.
!
!  N       (input) INTEGER
!          The order of the matrix T. N >= 0.
!
!  T       (input/output) REAL array, dimension (LDT,N)
!          On entry, the upper quasi-triangular matrix T, in Schur
!          canonical form.
!          On exit, T is overwritten by the reordered matrix T, again in
!          Schur canonical form, with the selected eigenvalues in the
!          leading diagonal blocks.
!
!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= max(1,N).
!
!  Q       (input/output) REAL array, dimension (LDQ,N)
!          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
!          On exit, if COMPQ = 'V', Q has been postmultiplied by the
!          orthogonal transformation matrix which reorders T; the
!          leading M columns of Q form an orthonormal basis for the
!          specified invariant subspace.
!          If COMPQ = 'N', Q is not referenced.
!
!  LDQ     (input) INTEGER
!          The leading dimension of the array Q.
!          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
!
!  WR      (output) REAL array, dimension (N)
!  WI      (output) REAL array, dimension (N)
!          The real and imaginary parts, respectively, of the reordered
!          eigenvalues of T. The eigenvalues are stored in the same
!          order as on the diagonal of T, with WR(i) = T(i,i) and, if
!          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
!          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
!          sufficiently ill-conditioned, then its value may differ
!          significantly from its value before reordering.
!
!  M       (output) INTEGER
!          The dimension of the specified invariant subspace.
!          0 < = M <= N.
!
!  S       (output) REAL
!          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
!          condition number for the selected cluster of eigenvalues.
!          S cannot underestimate the true reciprocal condition number
!          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
!          If JOB = 'N' or 'V', S is not referenced.
!
!  SEP     (output) REAL
!          If JOB = 'V' or 'B', SEP is the estimated reciprocal
!          condition number of the specified invariant subspace. If
!          M = 0 or N, SEP = norm(T).
!          If JOB = 'N' or 'E', SEP is not referenced.
!
!  WORK    (workspace) REAL array, dimension (LWORK)
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If JOB = 'N', LWORK >= max(1,N);
!          if JOB = 'E', LWORK >= M*(N-M);
!          if JOB = 'V' or 'B', LWORK >= 2*M*(N-M).
!
!  IWORK   (workspace) INTEGER array, dimension (LIWORK)
!          IF JOB = 'N' or 'E', IWORK is not referenced.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          If JOB = 'N' or 'E', LIWORK >= 1;
!          if JOB = 'V' or 'B', LIWORK >= M*(N-M).
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          = 1: reordering of T failed because some eigenvalues are too
!               close to separate (the problem is very ill-conditioned);
!               T may have been partially reordered, and WR and WI
!               contain the eigenvalues in the same order as in T; S and
!               SEP (if requested) are set to zero.
!
!  Further Details
!  ===============
!
!  STRSEN first collects the selected eigenvalues by computing an
!  orthogonal transformation Z to move them to the top left corner of T.
!  In other words, the selected eigenvalues are the eigenvalues of T11
!  in:
!
!                Z'*T*Z = ( T11 T12 ) n1
!                         (  0  T22 ) n2
!                            n1  n2
!
!  where N = n1+n2 and Z' means the transpose of Z. The first n1 columns
!  of Z span the specified invariant subspace of T.
!
!  If T has been obtained from the real Schur factorization of a matrix
!  A = Q*T*Q', then the reordered real Schur factorization of A is given
!  by A = (Q*Z)*(Z'*T*Z)*(Q*Z)', and the first n1 columns of Q*Z span
!  the corresponding invariant subspace of A.
!
!  The reciprocal condition number of the average of the eigenvalues of
!  T11 may be returned in S. S lies between 0 (very badly conditioned)
!  and 1 (very well conditioned). It is computed as follows. First we
!  compute R so that
!
!                         P = ( I  R ) n1
!                             ( 0  0 ) n2
!                               n1 n2
!
!  is the projector on the invariant subspace associated with T11.
!  R is the solution of the Sylvester equation:
!
!                        T11*R - R*T22 = T12.
!
!  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
!  the two-norm of M. Then S is computed as the lower bound
!
!                      (1 + F-norm(R)**2)**(-1/2)
!
!  on the reciprocal of 2-norm(P), the true reciprocal condition number.
!  S cannot underestimate 1 / 2-norm(P) by more than a factor of
!  sqrt(N).
!
!  An approximate error bound for the computed average of the
!  eigenvalues of T11 is
!
!                         EPS * norm(T) / S
!
!  where EPS is the machine precision.
!
!  The reciprocal condition number of the right invariant subspace
!  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
!  SEP is defined as the separation of T11 and T22:
!
!                     sep( T11, T22 ) = sigma-min( C )
!
!  where sigma-min(C) is the smallest singular value of the
!  n1*n2-by-n1*n2 matrix
!
!     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
!
!  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
!  product. We estimate sigma-min(C) by the reciprocal of an estimate of
!  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
!  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
!
!  When SEP is small, small changes in T can cause large changes in
!  the invariant subspace. An approximate bound on the maximum angular
!  error in the computed right invariant subspace is
!
!                      EPS * norm(T) / SEP
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            PAIR, SWAP, WANTBH, WANTQ, WANTS, WANTSP
      INTEGER            IERR, K, KASE, KK, KS, N1, N2, NN
      REAL               EST, RNORM, SCALE
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL               SLANGE
      EXTERNAL           LSAME, SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLACON, SLACPY, STREXC, STRSYL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      WANTBH = LSAME( JOB, 'B' )
      WANTS = LSAME( JOB, 'E' ) .OR. WANTBH
      WANTSP = LSAME( JOB, 'V' ) .OR. WANTBH
      WANTQ = LSAME( COMPQ, 'V' )
!
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.WANTS .AND. .NOT.WANTSP ) &
           THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.WANTQ ) THEN
         INFO = -2
      ELSE IF( N<0 ) THEN
         INFO = -4
      ELSE IF( LDT<MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDQ<1 .OR. ( WANTQ .AND. LDQ<N ) ) THEN
         INFO = -8
      ELSE
!
!        Set M to the dimension of the specified invariant subspace,
!        and test LWORK and LIWORK.
!
         M = 0
         PAIR = .FALSE.
         DO 10 K = 1, N
            IF( PAIR ) THEN
               PAIR = .FALSE.
            ELSE
               IF( K<N ) THEN
                  IF( T( K+1, K )==ZERO ) THEN
                     IF( SELECT( K ) ) &
                        M = M + 1
                  ELSE
                     PAIR = .TRUE.
                     IF( SELECT( K ) .OR. SELECT( K+1 ) ) &
                        M = M + 2
                  END IF
               ELSE
                  IF( SELECT( N ) ) &
                     M = M + 1
               END IF
            END IF
   10    CONTINUE
!
         N1 = M
         N2 = N - M
         NN = N1*N2
!
         IF( LWORK<1 .OR. ( ( WANTS .AND. .NOT.WANTSP ) .AND. &
             LWORK<NN ) .OR. ( WANTSP .AND. LWORK<2*NN ) ) THEN
            INFO = -15
         ELSE IF( LIWORK<1 .OR. ( WANTSP .AND. LIWORK<NN ) ) THEN
            INFO = -17
         END IF
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'STRSEN', -INFO )
         RETURN
      END IF
!
!     Quick return if possible.
!
      IF( M==N .OR. M==0 ) THEN
         IF( WANTS ) &
            S = ONE
         IF( WANTSP ) &
            SEP = SLANGE( '1', N, N, T, LDT, WORK )
         GO TO 40
      END IF
!
!     Collect the selected blocks at the top-left corner of T.
!
      KS = 0
      PAIR = .FALSE.
      DO 20 K = 1, N
         IF( PAIR ) THEN
            PAIR = .FALSE.
         ELSE
            SWAP = SELECT( K )
            IF( K<N ) THEN
               IF( T( K+1, K )/=ZERO ) THEN
                  PAIR = .TRUE.
                  SWAP = SWAP .OR. SELECT( K+1 )
               END IF
            END IF
            IF( SWAP ) THEN
               KS = KS + 1
!
!              Swap the K-th block to position KS.
!
               IERR = 0
               KK = K
               IF( K/=KS ) &
                  CALL STREXC( COMPQ, N, T, LDT, Q, LDQ, KK, KS, WORK, &
                               IERR )
               IF( IERR==1 .OR. IERR==2 ) THEN
!
!                 Blocks too close to swap: exit.
!
                  INFO = 1
                  IF( WANTS ) &
                     S = ZERO
                  IF( WANTSP ) &
                     SEP = ZERO
                  GO TO 40
               END IF
               IF( PAIR ) &
                  KS = KS + 1
            END IF
         END IF
   20 CONTINUE
!
      IF( WANTS ) THEN
!
!        Solve Sylvester equation for R:
!
!           T11*R - R*T22 = scale*T12
!
         CALL SLACPY( 'F', N1, N2, T( 1, N1+1 ), LDT, WORK, N1 )
         CALL STRSYL( 'N', 'N', -1, N1, N2, T, LDT, T( N1+1, N1+1 ), &
                      LDT, WORK, N1, SCALE, IERR )
!
!        Estimate the reciprocal of the condition number of the cluster
!        of eigenvalues.
!
         RNORM = SLANGE( 'F', N1, N2, WORK, N1, WORK )
         IF( RNORM==ZERO ) THEN
            S = ONE
         ELSE
            S = SCALE / ( SQRT( SCALE*SCALE / RNORM+RNORM )* &
                SQRT( RNORM ) )
         END IF
      END IF
!
      IF( WANTSP ) THEN
!
!        Estimate sep(T11,T22).
!
         EST = ZERO
         KASE = 0
   30    CONTINUE
         CALL SLACON( NN, WORK( NN+1 ), WORK, IWORK, EST, KASE )
         IF( KASE/=0 ) THEN
            IF( KASE==1 ) THEN
!
!              Solve  T11*R - R*T22 = scale*X.
!
               CALL STRSYL( 'N', 'N', -1, N1, N2, T, LDT, &
                            T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, &
                            IERR )
            ELSE
!
!              Solve  T11'*R - R*T22' = scale*X.
!
               CALL STRSYL( 'T', 'T', -1, N1, N2, T, LDT, &
                            T( N1+1, N1+1 ), LDT, WORK, N1, SCALE, &
                            IERR )
            END IF
            GO TO 30
         END IF
!
         SEP = SCALE / EST
      END IF
!
   40 CONTINUE
!
!     Store the output eigenvalues in WR and WI.
!
      DO 50 K = 1, N
         WR( K ) = T( K, K )
         WI( K ) = ZERO
   50 CONTINUE
      DO 60 K = 1, N - 1
         IF( T( K+1, K )/=ZERO ) THEN
            WI( K ) = SQRT( ABS( T( K, K+1 ) ) )* &
                      SQRT( ABS( T( K+1, K ) ) )
            WI( K+1 ) = -WI( K )
         END IF
   60 CONTINUE
      RETURN
!
!     End of STRSEN
!
      END
      SUBROUTINE STRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C, &
                         LDC, SCALE, INFO )
!
!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     March 31, 1993
!
!     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      REAL               SCALE
!     ..
!     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
!     ..
!
!  Purpose
!  =======
!
!  STRSYL solves the real Sylvester matrix equation:
!
!     op(A)*X + X*op(B) = scale*C or
!     op(A)*X - X*op(B) = scale*C,
!
!  where op(A) = A or A**T, and  A and B are both upper quasi-
!  triangular. A is M-by-M and B is N-by-N; the right hand side C and
!  the solution X are M-by-N; and scale is an output scale factor, set
!  <= 1 to avoid overflow in X.
!
!  A and B must be in Schur canonical form (as returned by SHSEQR), that
!  is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
!  each 2-by-2 diagonal block has its diagonal elements equal and its
!  off-diagonal elements of opposite sign.
!
!  Arguments
!  =========
!
!  TRANA   (input) CHARACTER*1
!          Specifies the option op(A):
!          = 'N': op(A) = A    (No transpose)
!          = 'T': op(A) = A**T (Transpose)
!          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
!
!  TRANB   (input) CHARACTER*1
!          Specifies the option op(B):
!          = 'N': op(B) = B    (No transpose)
!          = 'T': op(B) = B**T (Transpose)
!          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
!
!  ISGN    (input) INTEGER
!          Specifies the sign in the equation:
!          = +1: solve op(A)*X + X*op(B) = scale*C
!          = -1: solve op(A)*X - X*op(B) = scale*C
!
!  M       (input) INTEGER
!          The order of the matrix A, and the number of rows in the
!          matrices X and C. M >= 0.
!
!  N       (input) INTEGER
!          The order of the matrix B, and the number of columns in the
!          matrices X and C. N >= 0.
!
!  A       (input) REAL array, dimension (LDA,M)
!          The upper quasi-triangular matrix A, in Schur canonical form.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,M).
!
!  B       (input) REAL array, dimension (LDB,N)
!          The upper quasi-triangular matrix B, in Schur canonical form.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B. LDB >= max(1,N).
!
!  C       (input/output) REAL array, dimension (LDC,N)
!          On entry, the M-by-N right hand side matrix C.
!          On exit, C is overwritten by the solution matrix X.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M)
!
!  SCALE   (output) REAL
!          The scale factor, scale, set <= 1 to avoid overflow in X.
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!          = 1: A and B have common or very close eigenvalues; perturbed
!               values were used to solve the equation (but the matrices
!               A and B are unchanged).
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT
      REAL               A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN, &
                         SMLNUM, SUML, SUMR, XNORM
!     ..
!     .. Local Arrays ..
      REAL               DUM( 1 ), VEC( 2, 2 ), X( 2, 2 )
!     ..
!     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT, SLAMCH, SLANGE
      EXTERNAL           LSAME, SDOT, SLAMCH, SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL           SLABAD, SLALN2, SLASY2, SSCAL, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, REAL
!     ..
!     .. Executable Statements ..
!
!     Decode and Test input parameters
!
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
!
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT. &
          LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT. &
               LSAME( TRANB, 'C' ) ) THEN
         INFO = -2
      ELSE IF( ISGN/=1 .AND. ISGN/=-1 ) THEN
         INFO = -3
      ELSE IF( M<0 ) THEN
         INFO = -4
      ELSE IF( N<0 ) THEN
         INFO = -5
      ELSE IF( LDA<MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB<MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC<MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO/=0 ) THEN
         CALL XERBLA( 'STRSYL', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M==0 .OR. N==0 ) &
         RETURN
!
!     Set constants to control overflow
!
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*REAL( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
!
      SMIN = MAX( SMLNUM, EPS*SLANGE( 'M', M, M, A, LDA, DUM ), &
             EPS*SLANGE( 'M', N, N, B, LDB, DUM ) )
!
      SCALE = ONE
      SGN = ISGN
!
      IF( NOTRNA .AND. NOTRNB ) THEN
!
!        Solve    A*X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-left corner column by column by
!
!         A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                  M                         L-1
!        R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)].
!                I=K+1                       J=1
!
!        Start column loop (index = L)
!        L1 (L2) : column index of the first (first) row of X(K,L).
!
         LNEXT = 1
         DO 70 L = 1, N
            IF( L<LNEXT ) &
               GO TO 70
            IF( L==N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L )/=ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L).
!
            KNEXT = M
            DO 60 K = M, 1, -1
               IF( K>KNEXT ) &
                  GO TO 60
               IF( K==1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 )/=ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
!
               IF( L1==L2 .AND. K1==K2 ) THEN
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                               C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11<=SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11<ONE .AND. DB>ONE ) THEN
                     IF( DB>BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  IF( SCALOC/=ONE ) THEN
                     DO 10 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
!
               ELSE IF( L1==L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 20 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1==K2 ) THEN
!
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                               C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                               C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 40 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL SLASY2( .FALSE., .FALSE., ISGN, 2, 2, &
                               A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC, &
                               2, SCALOC, X, 2, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 50 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   50                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
!
   60       CONTINUE
!
   70    CONTINUE
!
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
!
!        Solve    A' *X + ISGN*X*B = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        upper-left corner column by column by
!
!          A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
!
!        Where
!                   K-1                        L-1
!          R(K,L) = SUM [A(I,K)'*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)]
!                   I=1                        J=1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         LNEXT = 1
         DO 130 L = 1, N
            IF( L<LNEXT ) &
               GO TO 130
            IF( L==N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L )/=ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
            KNEXT = 1
            DO 120 K = 1, M
               IF( K<KNEXT ) &
                  GO TO 120
               IF( K==M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K )/=ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
!
               IF( L1==L2 .AND. K1==K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11<=SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11<ONE .AND. DB>ONE ) THEN
                     IF( DB>BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  IF( SCALOC/=ONE ) THEN
                     DO 80 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
!
               ELSE IF( L1==L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 90 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1==K2 ) THEN
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 100 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL SLASY2( .TRUE., .FALSE., ISGN, 2, 2, A( K1, K1 ), &
                               LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, &
                               2, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 110 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  110                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
!
  120       CONTINUE
  130    CONTINUE
!
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
!
!        Solve    A'*X + ISGN*X*B' = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        top-right corner column by column by
!
!           A(K,K)'*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
!
!        Where
!                     K-1                          N
!            R(K,L) = SUM [A(I,K)'*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
!                     I=1                        J=L+1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         LNEXT = N
         DO 190 L = N, 1, -1
            IF( L>LNEXT ) &
               GO TO 190
            IF( L==1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 )/=ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
            KNEXT = 1
            DO 180 K = 1, M
               IF( K<KNEXT ) &
                  GO TO 180
               IF( K==M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K )/=ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
!
               IF( L1==L2 .AND. K1==K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC, &
                               B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11<=SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11<ONE .AND. DB>ONE ) THEN
                     IF( DB>BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  IF( SCALOC/=ONE ) THEN
                     DO 140 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  140                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
!
               ELSE IF( L1==L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 150 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  150                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1==K2 ) THEN
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 160 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  160                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                               B( L2, MIN(L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL SLASY2( .TRUE., .TRUE., ISGN, 2, 2, A( K1, K1 ), &
                               LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, &
                               2, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 170 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  170                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
!
  180       CONTINUE
  190    CONTINUE
!
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
!
!        Solve    A*X + ISGN*X*B' = scale*C.
!
!        The (K,L)th block of X is determined starting from
!        bottom-right corner column by column by
!
!            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L)' = C(K,L) - R(K,L)
!
!        Where
!                      M                          N
!            R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B(L,J)'].
!                    I=K+1                      J=L+1
!
!        Start column loop (index = L)
!        L1 (L2): column index of the first (last) row of X(K,L)
!
         LNEXT = N
         DO 250 L = N, 1, -1
            IF( L>LNEXT ) &
               GO TO 250
            IF( L==1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 )/=ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
!
!           Start row loop (index = K)
!           K1 (K2): row index of the first (last) row of X(K,L)
!
            KNEXT = M
            DO 240 K = M, 1, -1
               IF( K>KNEXT ) &
                  GO TO 240
               IF( K==1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 )/=ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
!
               IF( L1==L2 .AND. K1==K2 ) THEN
                  SUML = SDOT( M-K1, A( K1, MIN(K1+1, M ) ), LDA, &
                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC, &
                               B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
!
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11<=SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11<ONE .AND. DB>ONE ) THEN
                     IF( DB>BIGNUM*DA11 ) &
                        SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
!
                  IF( SCALOC/=ONE ) THEN
                     DO 200 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  200                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
!
               ELSE IF( L1==L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ), &
                               LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 210 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  210                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1==K2 ) THEN
!
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                               C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
!
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA, &
                               C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
!
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ), &
                               LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ), &
                               ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 220 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  220                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
!
               ELSE IF( L1/=L2 .AND. K1/=K2 ) THEN
!
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC, &
                               B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                               B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
!
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA, &
                               C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC, &
                               B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
!
                  CALL SLASY2( .FALSE., .TRUE., ISGN, 2, 2, A( K1, K1 ), &
                               LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X, &
                               2, XNORM, IERR )
                  IF( IERR/=0 ) &
                     INFO = 1
!
                  IF( SCALOC/=ONE ) THEN
                     DO 230 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  230                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
!
  240       CONTINUE
  250    CONTINUE
!
      END IF
!
      RETURN
!
!     End of STRSYL
!
      END
      SUBROUTINE XERBLA( SRNAME, INFO )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!     .. Scalar Arguments ..
      CHARACTER*6        SRNAME
      INTEGER            INFO
!     ..
!
!  Purpose
!  =======
!
!  XERBLA  is an error handler for the LAPACK routines.
!  It is called by an LAPACK routine if an input parameter has an
!  invalid value.  A message is printed and execution stops.
!
!  Installers may consider modifying the STOP statement in order to
!  call system-specific exception-handling facilities.
!
!  Arguments
!  =========
!
!  SRNAME  (input) CHARACTER*6
!          The name of the routine which called XERBLA.
!
!  INFO    (input) INTEGER
!          The position of the invalid parameter in the parameter list
!          of the calling routine.
!
! =====================================================================
!
!     .. Executable Statements ..
!
      WRITE( *, FMT = 9999 )SRNAME, INFO
!
      STOP
!
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ', &
            'an illegal value' )
!
!     End of XERBLA
!
      END
      SUBROUTINE XLAENV( ISPEC, NVALUE )
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992
!
!     .. Scalar Arguments ..
      INTEGER            ISPEC, NVALUE
!     ..
!
!  Purpose
!  =======
!
!  XLAENV sets certain machine- and problem-dependent quantities
!  which will later be retrieved by ILAENV.
!
!  Arguments
!  =========
!
!  ISPEC   (input) INTEGER
!          Specifies the parameter to be set in the COMMON array IPARMS.
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
!               the matrix to a triangular form)
!          = 7: the number of processors
!          = 8: another crossover point, for the multishift QR and QZ
!               methods for nonsymmetric eigenvalue problems.
!
!  NVALUE  (input) INTEGER
!          The value of the parameter specified by ISPEC.
!
!  =====================================================================
!
!     .. Arrays in Common ..
      INTEGER            IPARMS( 100 )
!     ..
!     .. Common blocks ..
      COMMON             / CLAENV / IPARMS
!     ..
!     .. Save statement ..
      SAVE               / CLAENV /
!     ..
!     .. Executable Statements ..
!
      IF( ISPEC>=1 .AND. ISPEC<=8 ) THEN
         IPARMS( ISPEC ) = NVALUE
      END IF
!
      RETURN
!
!     End of XLAENV
!
      END
