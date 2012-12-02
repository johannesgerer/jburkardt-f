! This file contains the test examples for the 3-dimensional cube
! that appear in the DECUHR paper:
! Espelid & Genz, Numerical Algorithms Vol 8, 1994.
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
REAL(kind=stnd), DIMENSION(1:NUMFUN) :: Value
REAL(kind=stnd) :: r
SELECT CASE (Fnr)
CASE(1)
 Value(1) = exp(x(1)+x(1)*x(2)+x(3)/3) /(sqrt(x(1)+x(2))*2.7878925361_stnd)
CASE(2)
 r = sqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))
 Value(1) = -log(r)*exp(x(1)*x(2)+x(3))/(sqrt(r)*0.1176364548_stnd)
END SELECT

RETURN
END FUNCTION F

SUBROUTINE INIT (invoer)
INTEGER, INTENT(IN) :: invoer
Fnr = invoer
RETURN
END SUBROUTINE INIT

END MODULE Integrand



PROGRAM Example3D

USE Precision_Model
USE CUI                      ! Cubpack User Interface
USE Integrand

INTEGER, PARAMETER :: n=3, & ! the dimension
                      m=1, & ! the number of simple regions
                      l=1, & ! the length of the integrand vector
                      Simplex = 1, &
                      Cube=2

INTEGER, DIMENSION(1:m) :: RgType
INTEGER :: NEval, i,j,k, testfunction, Job, Key
INTEGER, DIMENSION(1:3) :: jobtypes = (/ 1,12,2 /)
INTEGER, DIMENSION(1:3) :: keytypes = (/ 0,3,4 /)
REAL(kind=stnd), DIMENSION(1:n,0:n,1:m) :: Vertices
REAL(kind=stnd), DIMENSION(1:l) ::  IntegralValue, AbsErr
REAL(kind=stnd) ::  EpsRel
LOGICAL :: RESTART


RgType(1) = Cube
Vertices(1:n,0,1) = (/0 , 0 , 0 /)
Vertices(1:n,1,1) = (/1 , 0 , 0 /)
Vertices(1:n,2,1) = (/0 , 1 , 0 /)
Vertices(1:n,3,1) = (/0 , 0 , 1 /)

do testfunction = 1,2
   Print *,"Testfunction ",testfunction
   Print *,"==============="
   CALL INIT(testfunction)
   do i = 1,size(jobtypes)
      Job = jobtypes(i)
      do k = 1,size(keytypes)
         Key = keytypes(k)
         Print *,"JOB = ",Job," and Key=",Key
         Print *,"  Err Req     Integral          Error est    Cost"
         epsrel = 0.1_stnd
         RESTART= .false.
         do j = 1,4
            CALL CUBATR(n,l,F,m,Vertices(:,:,1:m),RgType,IntegralValue,&
              AbsErr,NEval=NEval,EpsRel=epsrel, &
              JOB=Job,KEY=Key,Restart=RESTART,MaxPTS=100000)
            Print '(E10.3,F22.15,E10.3,I8)',EpsRel, IntegralValue, AbsErr, NEval
            epsrel=epsrel*0.01_stnd
            if (Job /= 2) then
               RESTART = .true.
            end if
         end do
         CALL CUBATR()
         Print *,""
      end do
      Print *,"----------------------------------------------------------"
   end do
end do
STOP
END PROGRAM Example3D
