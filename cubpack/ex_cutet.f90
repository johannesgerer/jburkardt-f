! This file contains the first test program of DCUTET,
! a routine for integration over tetrahedrons.
! See Berntsen, Cools & Espelid, ACM TOMS Vol 19, 1993.
!
MODULE Integrand
IMPLICIT NONE
PUBLIC :: F
CONTAINS
FUNCTION F(NUMFUN,X) RESULT(Value)
USE Precision_Model
INTEGER, INTENT(IN) :: NUMFUN
REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
  Value(1) = EXP(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))
RETURN
END FUNCTION F

END MODULE Integrand



PROGRAM Example_CUTET

USE Precision_Model
USE CUI                      ! Cubpack User Interface
USE Integrand
IMPLICIT NONE

INTEGER, PARAMETER :: n=3, & ! the dimension
                      m=1, & ! the number of simple regions
                      l=1, & ! the length of the integrand vector
                      Tetrahedron = 1, &
                      Parallelogram = 2

INTEGER, DIMENSION(1:m) :: RgType
INTEGER :: NEval
REAL(kind=stnd), DIMENSION(1:n,0:n,1:m) :: Vertices
REAL(kind=stnd), DIMENSION(1:l) ::  IntegralValue, AbsErr
REAL(kind=stnd) :: epsrel


RgType(1) = Tetrahedron
Vertices(1:n,0,1) = (/0 , 0, 0 /)
Vertices(1:n,1,1) = (/1 , 0, 0 /)
Vertices(1:n,2,1) = (/0 , 1, 0 /)
Vertices(1:n,3,1) = (/0 , 0, 1 /)

WRITE(unit=*,fmt=*) "Simulation of DCUTET:"
epsrel = 1.0e-5_stnd
WRITE(unit=*,fmt=*) "    CUBATR will now be called with epsrel = ",epsrel

CALL CUBATR(n,l,F,m,Vertices,RgType,IntegralValue,AbsErr,Epsrel=epsrel,NEval=NEval,JOB=11)
WRITE(unit=*,fmt=*) "      Integral = ",IntegralValue
WRITE(unit=*,fmt=*) "      with estimated error < ",AbsErr
WRITE(unit=*,fmt=*) "      The number of integrand evaluations used = ",NEval

epsrel = 1.0e-8_stnd
WRITE(unit=*,fmt=*) "    CUBATR will now be called with epsrel = ",epsrel

CALL CUBATR(n,l,F,m,Vertices,RgType,IntegralValue,AbsErr,Epsrel=epsrel,NEval=NEval,Restart=.true.,JOB=11)
WRITE(unit=*,fmt=*) "      Integral = ",IntegralValue
WRITE(unit=*,fmt=*) "      with estimated error < ",AbsErr
WRITE(unit=*,fmt=*) "      The number of additional integrand evaluations used = ",NEval
STOP
END PROGRAM Example_CUTET
