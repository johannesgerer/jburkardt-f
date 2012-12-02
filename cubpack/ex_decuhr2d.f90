! This file contains a test example for the square
! that appears in the DECUHR paper:
! Espelid & Genz, Numerical Algorithms Vol 8, 1994.

! The restart feature is used with the default integration routine.

MODULE Integrand
IMPLICIT NONE
PUBLIC :: F
PRIVATE

CONTAINS
FUNCTION F(NUMFUN,X) RESULT(Value)
USE Precision_Model
INTEGER, INTENT(IN) :: NUMFUN
REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
REAL(kind=stnd), DIMENSION(1:NUMFUN) :: Value
! The integrand is scaled so that the exact result is 1.
 Value(1) = exp(2*x(1)+x(2)*(1-x(1))) * (1-x(1)) / sqrt(x(1))
 Value(1) = Value(1)/3.22289153891637_stnd
RETURN
END FUNCTION F

END MODULE Integrand



PROGRAM Example2D

USE Precision_Model
USE CUI                      ! Cubpack User Interface
USE Integrand

INTEGER, PARAMETER :: n=2, & ! the dimension
                      m=1, & ! the number of simple regions
                      l=1, & ! the length of the integrand vector
                      Cube=2

INTEGER, DIMENSION(1:m) :: RgType
INTEGER :: NEval, j
REAL(kind=stnd), DIMENSION(1:n,0:n,1:m) :: Vertices
REAL(kind=stnd), DIMENSION(1:l) ::  IntegralValue, AbsErr
REAL(kind=stnd) ::  EpsRel
LOGICAL :: Restart

RgType(1) = Cube
Vertices(1:n,0,1) = (/0 , 0  /)
Vertices(1:n,1,1) = (/1 , 0  /)
Vertices(1:n,2,1) = (/0 , 1  /)

 epsrel = 0.01_stnd
 Restart= .false.
 do j = 1,4
    Print *,"Request relative accuracy = ",epsrel
    CALL CUBATR(n,l,F,m,Vertices(:,:,1:m),RgType,IntegralValue,AbsErr, &
       NEval=NEval,EpsRel=epsrel,Restart=Restart)
    Print *,"-> Integral approximation = ",IntegralValue
    Print *,"-> with estimated error < ",AbsErr
    Print *,"-> The number of integrand evaluations used = ",NEval
    Print *,""
    epsrel=epsrel*0.01_stnd
    Restart = .true.
 end do

STOP
END PROGRAM Example2D
