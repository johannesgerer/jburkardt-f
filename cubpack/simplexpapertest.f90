! This file contains the full example of the paper
! A. Genz & R. Cools
! An adaptive numerical cubature algorithm for simplices
!
MODULE INTEGRAND
  USE PRECISION_MODEL
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: F
CONTAINS
  FUNCTION F( L, X ) RESULT(FUN)
    INTEGER,                       INTENT(IN) :: L
    REAL(kind=stnd), DIMENSION(:), INTENT(IN) :: X
    REAL(kind=stnd), DIMENSION(0:L-1)         :: FUN
    !
    REAL(kind=stnd)    :: S
    INTEGER            :: I
    INTEGER, PARAMETER :: N = 5
    S = 0
    DO I = 1, N
       S = S + ( I*X(I) )**2
    END DO
    FUN(0) = EXP(-S)
    FUN(1:N) = X(1:N)*FUN(0)
  END FUNCTION F
!
END MODULE INTEGRAND
!
PROGRAM Simplex_Example

USE Precision_Model
USE CUI                                 ! Cubpack User Interface
USE INTEGRAND

INTEGER, PARAMETER         :: n = 5,   & ! the dimension
                              nf = n + 1 ! number of integrand functions

INTEGER                             :: NEval, i, Inform
REAL(kind=stnd), DIMENSION(n,0:n,1) :: Simplex
REAL(kind=stnd), DIMENSION(n,0:n,2) :: TwoSimplices
REAL(kind=stnd), DIMENSION(0:n)     :: Value, AbsErr

Inform = -1                      ! Soft noisy errors            

Simplex = 0                      ! Build the unit simplex 
DO i = 1, n
   Simplex(i,i,1) = 1
END DO
TwoSimplices(:,:,1:1) = Simplex    ! Split it into 2 parts.
TwoSimplices(1,0,1) = 0.5_stnd     !  Therefore one has to change one
TwoSimplices(:,:,2:2) = Simplex    !    coordinate of the unit simplex.
TwoSimplices(1,1,2) = 0.5_stnd

!   The two subregions | are simplices    | |
!                      v                  v v
CALL CUBATR( n, nf, F, 2, TwoSimplices, (/1,1/), Value, AbsErr,    &
             Inform, NEval, EpsRel=5.0e-4_stnd, MaxPts=200000 )

Print "(""Expected values are ""    /5F12.8)",  Value(1:n)/Value(0) 
Print "(""with estimated errors < ""/5F12.8)", AbsErr(1:n)/Value(0) 
Print *,"The number of integrand evaluations used was ", NEval 

Inform = -1                      ! Soft noisy errors            
!    The single region | is a simplex    |
!                      v                 v
CALL CUBATR( n, nf, F, 1, Simplex,     (/1/), Value, AbsErr,      &
             Inform, NEval, EpsRel=5.0e-4_stnd, MaxPts=200000 )

Print "(""Expected values are ""    /5F12.8)",  Value(1:n)/Value(0) 
Print "(""with estimated errors < ""/5F12.8)", AbsErr(1:n)/Value(0) 
Print *,"The number of integrand evaluations used was ", NEval 

END PROGRAM Simplex_Example
