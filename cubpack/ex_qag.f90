MODULE Integrand
PRIVATE
PUBLIC :: F
CONTAINS
   FUNCTION F(NUMFUN,X) RESULT(Value)
   USE Precision_Model
   INTEGER, INTENT(IN) :: NUMFUN
   REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
   REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
   Value(1) = (abs(x(1)-1.0_stnd/3.0_stnd))**(0.8_stnd)
   
   RETURN 
   END FUNCTION F
END MODULE Integrand

PROGRAM Example_QAG
! 
! This is equivalent to the test program for QAG from QUADPACK
! The output is identical to the orginal Fortran 77 programs.
! 25 May 1999

USE Precision_Model
USE CUI                      ! Cubpack User Interface
USE Integrand

IMPLICIT NONE

INTEGER, PARAMETER :: n=1, & ! the dimension
                      Finite_interval = 1

INTEGER ::   RgType, NEval,  j, key
REAL(kind=stnd), DIMENSION(1:n,0:n) :: Vertices
REAL(kind=stnd) :: IntegralValue, AbsErr, epsrel

epsrel = 0.1_stnd

do j=1,10
   epsrel = epsrel*0.1_stnd
   do key = 1,6
      RgType = Finite_interval
      Vertices(1,:) = (/ 0 , 1 /)
      CALL CUBATR(n,F,Vertices,RgType,IntegralValue,AbsErr, &
                  MAXPTS=50000,EpsRel=epsrel,Key=key,NEval=NEval,JOB=1)
      write(unit=*,fmt="( "" Results for key = "",I3 )") KEY
      write(unit=*,fmt="( "" Results for epsrel = "",es9.2 )") EPSREL
      write(unit=*,fmt="( "" INTEGRAL APPROXIMATION = "",es15.8 )") IntegralValue
      write(unit=*,fmt="( "" ESTIMATE OF ABSOLUTE ERROR = "",es9.2 )") ABSERR
      write(unit=*,fmt="( "" NUMBER OF FUNCTION EVALATIONS = "",i5 )") NEVAL

   end do
end do
STOP 
END PROGRAM Example_QAG                                          
