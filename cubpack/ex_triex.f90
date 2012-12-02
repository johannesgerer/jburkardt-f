!
! This file is equivalent to the tests mentioned in the TRIEX paper.
! The results are NOT fully identical because CUBPACK uses another
! integration rule and error estimator for a triangle.
!
MODULE Integrand
IMPLICIT NONE
PUBLIC :: F, INIT
PRIVATE
INTEGER, PRIVATE :: Fnr
CONTAINS
FUNCTION F(NUMFUN,X) RESULT(Value)
USE Precision_Model
INTEGER, INTENT(IN) :: NUMFUN
REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
REAL(kind=stnd), DIMENSION(NUMFUN) :: Value
INTEGER, DIMENSION(1:10), PARAMETER ::  &
                     P = (/ 1,0,1,1,0,0,1,1,1,1 /),         &
                     Q = (/ 0,0,0,0,0,0,1,1,1,0 /),         &
                     R = (/ 0,1,0,0,1,1,0,1,0,1 /)
REAL(kind=stnd), DIMENSION(1:10), PARAMETER :: &
         A=(/0.0_stnd,0.0_stnd,0.5_stnd,1.0/3.0_stnd,0.0_stnd,0.0_stnd,1.0_stnd,0.0_stnd,0.5_stnd,0.0_stnd /), &
         B = 0,                                               &
         C = (/0.0_stnd,1.0_stnd,0.0_stnd,0.0_stnd,0.5_stnd,1.0/3.0_stnd,0.0_stnd,0.0_stnd,0.0_stnd,1.0/3.0_stnd/)
SELECT CASE (Fnr)
CASE (1:10)
   Value(1) = p(Fnr)*sqrt(abs(x(1)-a(Fnr)))+q(Fnr)*sqrt(abs(x(2)-b(Fnr))) &
              + r(Fnr)*sqrt(abs(x(1)+x(2)-c(Fnr)))
CASE(11)
   Value(1) = log(x(1)+x(2))
CASE (12)
   Value(1) = 1/sqrt(x(1)*x(1)+x(2)*x(2))
CASE (13)
   Value(1) = log(sqrt(x(1)*x(1)+x(2)*x(2)))/sqrt(x(1)*x(1)+x(2)*x(2))
CASE (14)
   Value(1) = sin(x(1)) * cos(5*x(2))
CASE (15)
   ! I guess there is a typo in the triex paper and t must be 1.
   Value(1) = sin(11*x(1)) * cos(x(2))   
CASE (16)
   Value(1) = x(1)**(-0.2_stnd)*9.0_stnd/6.25_stnd
CASE (17)
   Value(1) = 1.0_stnd/sqrt(x(1)) + 1.0_stnd/sqrt(x(2))     &
             + 1.0_stnd/sqrt(x(1) + x(2))
   Value(1) = Value(1)*3.0_stnd/10.0_stnd
CASE (18)
   Value(1) = 3.0_stnd/(sqrt(1 - x(1) - x(2))*4.0_stnd)
CASE (19)
   Value(1) = (x(1)*x(2))**(-0.2_stnd)/0.9481026454955768_stnd
CASE (20)
   Value(1) = -2.0_stnd * log(x(1)*x(2))/ 3.0_stnd
CASE (21)
   Value(1) = (1.0_stnd/sqrt(abs(x(1)-0.25_stnd)) +          &
              1.0_stnd/sqrt(abs(x(2) - 0.5_stnd)))/3.11357229949_stnd
END SELECT
RETURN
END FUNCTION F

SUBROUTINE INIT (invoer)
INTEGER, INTENT(IN) :: invoer
Fnr = invoer
RETURN
END SUBROUTINE INIT

END MODULE Integrand

! ---- Now follows the main program ---

PROGRAM Example_Triex

USE Precision_Model
USE CUI                      ! Cubpack User Interface
USE Integrand
IMPLICIT NONE

INTEGER, PARAMETER :: n=2, & ! the dimension
                      m=1, & ! the number of simple regions
                      l=1, & ! the length of the integrand vector
                      Triangle = 1, &  ! increase readability
                      Parallelogram = 2

INTEGER, DIMENSION(1:n) :: RgType
INTEGER :: NEval,i
REAL(kind=stnd), DIMENSION(1:n,0:n,1:m) :: Vertices
REAL(kind=stnd), DIMENSION(1:l) ::  IntegralValue, AbsErr
REAL(kind=stnd) :: epsrel

RgType(1) = Triangle
Vertices(1:n,0,1) = (/0 , 0 /)
Vertices(1:n,1,1) = (/1 , 0 /)
Vertices(1:n,2,1) = (/0 , 1 /)

epsrel = 1.0e-6_stnd
do i = 1,21
   WRITE(unit=*,fmt=*) 
   WRITE(unit=*,fmt=*) "Testfunction ",i
   CALL INIT(i)

   WRITE(unit=*,fmt=*) "CUBATR will now be called with epsrel = ",real(epsrel)
   CALL CUBATR(n,l,F,m,Vertices(:,:,1:m),RgType,IntegralValue, &
        AbsErr,Epsrel=epsrel,NEval=NEval,JOB=2,MaxPts=100000)
   WRITE(unit=*,fmt=*) "Integral = ",IntegralValue
   WRITE(unit=*,fmt=*) "with estimated error < ",real(AbsErr)
   WRITE(unit=*,fmt=*) "The number of integrand evaluations used = ",NEval
   CALL CUBATR()
end do
STOP
END PROGRAM Example_Triex
