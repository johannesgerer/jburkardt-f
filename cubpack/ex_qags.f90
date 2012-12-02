MODULE Integrand
PRIVATE
PUBLIC :: F
CONTAINS
   FUNCTION F(NUMFUN,X) RESULT(Value)
   USE Precision_Model
   INTEGER, INTENT(IN) :: NUMFUN
   REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
   REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
   Value(1) = log(x(1))*(x(1)**0.2_stnd)
   RETURN 
   END FUNCTION F
END MODULE Integrand

PROGRAM Example_QAGS
! 
! This is equivalent to the test program for QAGS from QUADPACK
! The output is identical to the orginal Fortran 77 programs
! (at least on Linux using NAGWare f95 and f77).
! Minor differences for small numbers possible on e.g. Sun.

USE Precision_Model
USE CUI                      ! Cubpack User Interface
USE Integrand

INTEGER, PARAMETER :: n=1, & ! the dimension
                      Finite_interval = 1

INTEGER  ::   RgType, NEval, j, key
REAL(kind=stnd):: IntegralValue, AbsErr, epsrel
REAL(kind=stnd), DIMENSION(1:n,0:n) :: Vertices


epsrel = 0.1_stnd

do j=1,10
   epsrel = epsrel*0.1_stnd
   RgType = Finite_interval
   Vertices(1,:) = (/ 0 , 1 /)
   CALL CUBATR(n,F,Vertices,RgType,IntegralValue,AbsErr, &
                  EpsRel=epsrel,NEval=NEval,JOB=2)
!                 EpsRel=epsrel,Key=2,NEval=NEval,JOB=2)
!                               ^^/^^             ^^^/^
! QAGS uses dqk21 only  ---------/                  / 
!           (Key=2 is therefore made default.)     /
!      uses extrapolation ------------------------/
!
!  write(unit=*,fmt="( "" Results for key = "",I3 )") KEY
   write(unit=*,fmt="( "" Resuls for epsrel = "",es9.2 )") EPSREL
   write(unit=*,fmt="( "" INTEGRAL APPROXIMATION = "",es15.8 )") IntegralValue
   write(unit=*,fmt="( "" ESTIMATE OF ABSOLUTE ERROR = "",es9.2 )") ABSERR
   write(unit=*,fmt="( "" NUMBER OF FUNCTION EVALATIONS = "",i5 )") NEVAL
   write(unit=*,fmt="( "" ERROR CODE = 0"" )") ! IER

end do
STOP 
END PROGRAM Example_QAGS
