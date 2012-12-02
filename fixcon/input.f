*DECK AVNTST
      SUBROUTINE AVNTST (LUN, KPRINT, IPASS)
!***BEGIN PROLOGUE  AVNTST
!***PURPOSE  Quick check for AVINT.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (AVNTST-S, DAVNTS-D)
!***AUTHOR  (UNKNOWN)
!***ROUTINES CALLED  AVINT, R1MACH, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901205  Changed usage of R1MACH(3) to R1MACH(4).  (RWC)
!   910501  Added PURPOSE and TYPE records.  (WRB)
!   910708  Minor modifications in use of KPRINT.  (WRB)
!   920210  Code restructured and revised to test error returns for all
!           values of KPRINT.  (WRB)
!***END PROLOGUE  AVNTST
      DIMENSION X(501), Y(501)
      LOGICAL FATAL
!***FIRST EXECUTABLE STATEMENT  AVNTST
      if ( KPRINT >= 2) write (LUN,9000)
      IPASS = 1
      TOL = MAX(.0001E0,SQRT(R1MACH(4)))
      TOL1 = 1.0E-2*TOL
!
!     Perform first accuracy test.
!
      A = 0.0E0
      B = 5.0E0
      XINT = EXP(5.0D0) - 1.0D0
      N = 500
      RN1 = N - 1
      SQB = SQRT(B)
      DEL = 0.4E0*(B-A)/(N-1)
      DO I = 1,N
        X(I) = SQB*SQRT(A+(I-1)*(B-A)/RN1) + DEL
        Y(I) = EXP(X(I))
      end do
      call AVINT (X, Y, N, A, B, ANS, IERR)
!
!  See if test was passed.
!
      if ( ABS(ANS-XINT)  >  TOL ) then
        IPASS = 0
        if ( KPRINT >= 3) write (LUN,9010) IERR, ANS, XINT
      end if
!
!  Perform second accuracy test.
!
      X(1) = 0.0E0
      X(2) = 5.0E0
      Y(1) = 1.0E0
      Y(2) = 0.5E0
      A = -0.5E0
      B = 0.5E0
      XINT = 1.0E0
      call AVINT (X, Y, 2, A, B, ANS, IERR)
!
!  See if test was passed.
!
      if ( ABS(ANS-XINT)  >  TOL1 ) then
        IPASS = 0
        if ( KPRINT >= 3) write (LUN,9010) IERR, ANS, XINT
      end if
!
!  Send message indicating passage or failure of tests.
!
      if ( KPRINT >= 2 ) then
        if ( IPASS  ==  1 ) then
           if ( KPRINT >= 3) write (LUN,9020)
        else
           write (LUN,9030)
        end if
      end if
!
!  Test error returns.
!
      call XGETF (KONTRL)
      if ( KPRINT  <=  2 ) then
         call XSETF (0)
      else
         call XSETF (1)
      end if
      FATAL = .FALSE.
      call XERCLR

      if ( KPRINT >= 3 ) then
        write (LUN,9040)
      end if
      DO 110 I = 1,20
        X(I) = (I-1)/19.0E0 - 0.01E0
        if ( I  /=  1) Y(I) = X(I)/(EXP(X(I))-1.0)
  110 CONTINUE
!
!  Test IERR = 1 error return.
!
      Y(1) = 1.0E0
      call AVINT (X, Y, 20, 0.0E0, 1.0E0, ANS, IERR)
      if ( IERR  /=  1 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9060) IERR, 1
      end if
      call XERCLR
!
!  Test IERR = 2 error return.
!
      call AVINT (X, Y, 20, 1.0E0, 0.0E0, ANS, IERR)
      if ( IERR  /=  2 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9060) IERR, 2
      end if
      if ( ANS  /=  0.0E0 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9070)
      end if
      call XERCLR
!
!  Test IERR = 5 error return.
!
      call AVINT (X, Y, 1, 0.0E0, 1.0E0, ANS, IERR)
      if ( IERR  /=  5 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9060) IERR, 5
      end if
      if ( ANS  /=  0.0E0 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9070)
      end if
      call XERCLR
!
!  Test IERR = 4 error return.
!
      X(1) = 1.0E0/19.0E0
      X(2) = 0.0E0
      call AVINT (X, Y, 20, 0.0E0, 1.0E0, ANS, IERR)
      if ( IERR  /=  4 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9060) IERR, 4
      end if
      if ( ANS  /=  0.0E0 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9070)
      end if
      call XERCLR
!
!  Test IERR = 3 error return.
!
      X(1) = 0.0E0
      X(2) = 1.0E0/19.0E0
      call AVINT (X, Y, 20, 0.0E0, .01E0, ANS, IERR)
      if ( IERR  /=  3 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9060) IERR, 3
      end if
      if ( ANS  /=  0.0E0 ) then
        IPASS = 0
        FATAL = .TRUE.
        if ( KPRINT >= 3) write (LUN,9070)
      end if
      call XERCLR
!
!  Reset XERMSG control variables and write summary.
!
      call XSETF (KONTRL)
      if ( FATAL ) then
         if ( KPRINT >= 2 ) then
            write (LUN, 9080)
         end if
      else
         if ( KPRINT >= 3 ) then
            write (LUN, 9090)
         end if
      end if
!
!  Write PASS/FAIL message.
!
      if ( IPASS == 1  .and.  KPRINT>=3) write (LUN,9100)
      if ( IPASS == 0  .and.  KPRINT>=2) write (LUN,9110)
      RETURN
 9000 FORMAT ('1' / ' AVINT Quick Check')
 9010 FORMAT (/' FAILED ACCURACY TEST' /
     +        ' IERR=', I2, 5X, 'COMPUTED ANS=', E20.11 / 14X,
     +        'CORRECT ANS=', E20.11, 5X, 'REQUESTED ERR=', E10.2)
 9020 FORMAT (/ ' AVINT passed both accuracy tests.')
 9030 FORMAT (/ ' AVINT failed at least one accuracy test.')
 9040 FORMAT (/ ' Test error returns from AVINT' /
     +        ' 4 error messages expected' /)
 9060 FORMAT (/' IERR =', I2, ' and it should =', I2 /)
 9070 FORMAT (1X, 'ANS  /=  0')
 9080 FORMAT (/ ' At least one incorrect argument test FAILED')
 9090 FORMAT (/ ' All incorrect argument tests PASSED')
 9100 FORMAT (/' ***************AVINT PASSED ALL TESTS***************')
 9110 FORMAT (/' ***************AVINT FAILED SOME TESTS**************')
      END
*DECK BIKCK
      SUBROUTINE BIKCK (LUN, KPRINT, IPASS)
!***BEGIN PROLOGUE  BIKCK
!***PURPOSE  Quick check for BESI and BESK.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BIKCK-S, DBIKCK-D)
!***KEYWORDS  QUICK CHECK
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!   BIKCK is a quick check routine for BESI and BESK.  The main loops
!   evaluate the Wronskian and test the error.  Underflow and overflow
!   diagnostics are checked in addition to illegal arguments.
!
!***ROUTINES CALLED  BESI, BESK, NUMXER, R1MACH, XERCLR, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891004  Removed unreachable code.  (WRB)
!   891004  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   901013  Editorial changes, some restructing and modifications to
!           obtain more information when there is failure of the
!           Wronskian.  (RWC)
!   910501  Added PURPOSE and TYPE records.  (WRB)
!   910708  Code revised to test error returns for all values of
!           KPRINT.  (WRB)
!***END PROLOGUE  BIKCK
      INTEGER I, IX, K, KONTRL, KODE, LUN, M, N, NERR, NU, NW, NY
      REAL ALP, DEL, ER, FNU, FNUP, RX, TOL, X
      REAL FN(3), W(5), XX(5), Y(5)
      REAL R1MACH
      LOGICAL FATAL
!***FIRST EXECUTABLE STATEMENT  BIKCK
      if ( KPRINT >= 2) write (LUN,90000)
!
      IPASS = 1
      XX(1) = 0.49E0
      XX(2) = 1.3E0
      XX(3) = 5.3E0
      XX(4) = 13.3E0
      XX(5) = 21.3E0
      FN(1) = 0.095E0
      FN(2) = 0.70E0
      FN(3) = 0.0E0
      TOL = 500.0E0*MAX(R1MACH(4), 7.1E-15)
      DO 60 KODE=1,2
         DO 50 M=1,3
            DO 40 N=1,4
               DO 30 NU=1,4
                  FNU = FN(M) + 12*(NU-1)
                  DO 20 IX=1,5
                     if ( IX < 2  .and.  NU > 3) GO TO 20
                     X = XX(IX)
                     RX = 1.0E0/X
                     call BESI(X, FNU, KODE, N, Y, NY)
                     if ( NY /= 0) GO TO 20
                     call BESK(X, FNU, KODE, N, W, NW)
                     if ( NW /= 0) GO TO 20
                     FNUP = FNU + N
                     call BESI(X,FNUP,KODE,1,Y(N+1),NY)
                     if ( NY /= 0) GO TO 20
                     call BESK(X,FNUP,KODE,1,W(N+1),NW)
                     if ( NW /= 0) GO TO 20
                     DO 10 I=1,N
                        ER = Y(I+1)*W(I) + W(I+1)*Y(I) - RX
                        ER = ABS(ER)*X
                        if ( ER > TOL ) then
                           IPASS = 0
                           if ( KPRINT>=2) write (LUN,90010) KODE,M,N,
     *                        NU,IX,I,X,ER,TOL,
     *                        Y(I),Y(I+1),W(I),W(I+1)
                        end if
   10                CONTINUE
   20             CONTINUE
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
!
!  Check small values of X and order
!
      N = 2
      FNU = 1.0E0
      X = R1MACH(4)/100.0E0
      DO 80 I=1,3
         DO 70 KODE=1,2
            call BESI(X, FNU, KODE, N, Y, NY)
            call BESK(X, FNU, KODE, N, W, NW)
            ER = Y(2)*W(1) + W(2)*Y(1) - 1.0E0/X
            ER = ABS(ER)*X
            if ( ER > TOL ) then
               IPASS = 0
               if ( KPRINT>=2) write (LUN,90020) I,KODE,FNU,X,ER,TOL,
     +            Y(1),Y(2),W(1),W(2)
               GO TO 700
            end if
   70    CONTINUE

  700    FNU = R1MACH(4)/100.0E0
         X = XX(2*I-1)
   80 CONTINUE
!
!  Check large values of X and order
!
      KODE = 2
      DO 76 K=1,2
         DEL = 30*(K-1)
         FNU = 45.0E0+DEL
         DO 75 N=1,2
            X = 20.0E0 + DEL
            DO 71 I=1,5
               RX = 1.0E0/X
               call BESI(X, FNU, KODE, N, Y, NY)
               if ( NY /= 0) GO TO 71
               call BESK(X, FNU, KODE, N, W, NW)
               if ( NW /= 0) GO TO 71
               if ( N == 1 ) then
                  FNUP = FNU + 1.0E0
                  call BESI(X,FNUP,KODE,1,Y(2),NY)
                  if ( NY /= 0) GO TO 71
                  call BESK(X,FNUP,KODE,1,W(2),NW)
                  if ( NW /= 0) GO TO 71
               end if
               ER = Y(2)*W(1) + Y(1)*W(2) - RX
               ER = ABS(ER)*X
               if ( ER > TOL ) then
                  IPASS = 0
                  if ( KPRINT>=2) write (LUN,90030) K,N,I,FNUP,X,
     +               ER,TOL,Y(1),Y(2),W(1),W(2)
                  GO TO 760
               end if
               X = X + 10.0E0
   71       CONTINUE
   75    CONTINUE
   76 CONTINUE
!
!  Check underflow flags
!
  760 X = R1MACH(1)*10.0E0
      ALP = 12.3E0
      N = 3
      call BESI(X, ALP, 1, N, Y, NY)
      if ( NY /= 3 ) then
         IPASS = 0
         if ( KPRINT>=2) write (LUN,90040)
      end if

      X = LOG(R1MACH(2)/10.0E0) + 20.0E0
      ALP = 1.3E0
      N = 3
      call BESK(X, ALP, 1, N, W, NW)
      if ( NW /= 3 ) then
         IPASS = 0
         if ( KPRINT>=2) write (LUN,90050)
      end if
!
!  Trigger 10 error conditions
!
      call XGETF (KONTRL)
      if ( KPRINT  <=  2 ) then
         call XSETF (0)
      else
         call XSETF (1)
      end if
      FATAL = .FALSE.
      call XERCLR

      if ( KPRINT >= 3) write (LUN,90060)
      XX(1) = 1.0E0
      XX(2) = 1.0E0
      XX(3) = 1.0E0
      XX(4) = 1.0E0
!
!  Illegal arguments
!
      DO I=1,4
         XX(I) = -XX(I)
         K = INT(XX(3))
         N = INT(XX(4))
         call BESI(XX(1), XX(2), K, N, Y, NY)
         if ( NUMXER(NERR)  /=  2 ) then
            IPASS = 0
            FATAL = .TRUE.
         end if
         call XERCLR
         call BESK(XX(1), XX(2), K, N, W, NW)
         if ( NUMXER(NERR)  /=  2 ) then
            IPASS = 0
            FATAL = .TRUE.
         end if
         call XERCLR
         XX(I) = -XX(I)
      end do
!
!     Trigger overflow
!
      X = LOG(R1MACH(2)/10.0E0) + 20.0E0
      N = 3
      ALP = 2.3E0
      call BESI(X, ALP, 1, N, Y, NY)
      if ( NUMXER(NERR)  /=  6 ) then
         IPASS = 0
         FATAL = .TRUE.
      end if
      call XERCLR

      X = R1MACH(1)*10.0E0
      call BESK(X, ALP, 1, N, W, NW)
      if ( NUMXER(NERR)  /=  6 ) then
         IPASS = 0
         FATAL = .TRUE.
      end if
      call XERCLR

      call XSETF (KONTRL)
      if ( FATAL ) then
         if ( KPRINT >= 2 ) then
            write (LUN, 90070)
         end if
      else
         if ( KPRINT >= 3 ) then
            write (LUN, 90080)
         end if
      end if

      if ( IPASS == 1  .and.  KPRINT>=2) write (LUN,90100)
      if ( IPASS == 0  .and.  KPRINT>=1) write (LUN,90110)
      RETURN

90000 FORMAT (/ ' QUICK CHECKS FOR BESI AND BESK' //)
90010 FORMAT (/ ' ERROR IN QUICK CHECK OF WRONSKIAN', 1P /
     +        ' KODE = ', I1,', M = ', I1, ', N = ', I1, ', NU = ', I1,
     +        ', IX = ', I1, ', I = ', I1 /
     +        ' X = ', E14.7, ', ER   = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(I) = ', E14.7, ', Y(I+1) = ', E14.7 /
     +        ' W(I) = ', E14.7, ', W(I+1) = ', E14.7)
90020 FORMAT (/ ' ERROR IN QUICK CHECK OF SMALL X AND ORDER', 1P /
     +        ' I = ', I1,', KODE = ', I1, ', FNU = ', E14.7 /
     +        ' X = ', E14.7, ', ER = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(1) = ', E14.7, ', Y(2) = ', E14.7 /
     +        ' W(1) = ', E14.7, ', W(2) = ', E14.7)
90030 FORMAT (/ ' ERROR IN QUICK CHECK OF LARGE X AND ORDER', 1P /
     +        ' K = ', I1,', N = ', I1, ', I = ', I1,
     +        ', FNUP = ', E14.7 /
     +        ' X = ', E14.7, ', ER = ', E14.7, ', TOL = ', E14.7 /
     +        ' Y(1) = ', E14.7, ', Y(2) = ', E14.7 /
     +        ' W(1) = ', E14.7, ', W(2) = ', E14.7)
90040 FORMAT (/ ' ERROR IN BESI UNDERFLOW TEST' /)
90050 FORMAT (/ ' ERROR IN BESK UNDERFLOW TEST' /)
90060 FORMAT (// ' TRIGGER 10 ERROR CONDITIONS' //)
90070 FORMAT (/ ' AT LEAST ONE INCORRECT ARGUMENT TEST FAILED')
90080 FORMAT (/ ' ALL INCORRECT ARGUMENT TESTS PASSED')
90100 FORMAT (/' **********BESI AND BESK PASSED ALL TESTS************')
90110 FORMAT (/' **********BESI OR BESK FAILED SOME TESTS************')
      END
