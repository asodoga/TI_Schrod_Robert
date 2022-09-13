SUBROUTINE Tana_F2_F1_Vep(F2,F1,Vep,Q)
USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64
IMPLICIT NONE


real(kind=real64), intent(in)    :: Q(3)
real(kind=real64), intent(inout) :: F1(3)
real(kind=real64), intent(inout) :: F2(3,3)
real(kind=real64), intent(inout) :: Vep

integer :: i,j



F2 = 0._real64
F1 = 0._real64

F2(2,3) = F2(2,3)  +0.00000784383625_real64 * Q(1)**(-1) * sin(Q(3))
F1(2)   = F1(2)    +0.00001176575438_real64 * Q(1)**(-1) * cos(Q(3))
F2(3,3) = F2(3,3)  +0.00001568767251_real64 * Q(1)**(-1) * Q(2)**(-1) * cos(Q(3))
F2(1,2) = F2(1,2)  -0.00001568767251_real64 * cos(Q(3))
F1(1)   = F1(1)    +0.00001176575438_real64 * Q(2)**(-1) * cos(Q(3))
F2(1,3) = F2(1,3)  +0.00000784383625_real64 * Q(2)**(-1) * sin(Q(3))**(-1)
F1(1)   = F1(1)    -0.00000392191813_real64 * Q(2)**(-1) * cos(Q(3))*sin(Q(3))**(-2)
F1(3)   = F1(3)    -0.00001568767251_real64 * Q(1)**(-1) * Q(2)**(-1) * sin(Q(3))**(-1)
F2(1,3) = F2(1,3)  -0.00000784383625_real64 * Q(2)**(-1) * cos(Q(3))**(2)*sin(Q(3))**(-1)
F1(1)   = F1(1)    +0.00000392191813_real64 * Q(2)**(-1) * cos(Q(3))**(3)*sin(Q(3))**(-2)
F1(3)   = F1(3)    +0.00001568767251_real64 * Q(1)**(-1) * Q(2)**(-1) * cos(Q(3))**(2)*sin(Q(3))**(-1)
F2(2,3) = F2(2,3)  +0.00000784383625_real64 * Q(1)**(-1) * sin(Q(3))**(-1)
F1(2)   = F1(2)    -0.00000392191813_real64 * Q(1)**(-1) * cos(Q(3))*sin(Q(3))**(-2)
F2(2,3) = F2(2,3)  -0.00000784383625_real64 * Q(1)**(-1) * cos(Q(3))**(2)*sin(Q(3))**(-1)
F1(2)   = F1(2)    +0.00000392191813_real64 * Q(1)**(-1) * cos(Q(3))**(3)*sin(Q(3))**(-2)
F2(1,3) = F2(1,3)  +0.00000784383625_real64 * Q(2)**(-1) * sin(Q(3))
F2(3,3) = F2(3,3)  -0.00028000412764_real64 * Q(1)**(-2)
F2(1,1) = F2(1,1)  -0.00028000412764_real64
F2(3,3) = F2(3,3)  -0.00014402858988_real64 * Q(2)**(-2)
F2(2,2) = F2(2,2)  -0.00014402858988_real64

DO i=1,3
DO j=i+1,3
  F2(j,i) = F2(i,j)
END DO
END DO

Vep = 0._real64
Vep = Vep -0.00000784383625_real64 * Q(1)**(-1) * Q(2)**(-1) * cos(Q(3))
Vep = Vep -0.00000392191813_real64 * Q(1)**(-1) * Q(2)**(-1) * cos(Q(3))**(3)*sin(Q(3))**(-2)
Vep = Vep +0.00000784383625_real64 * Q(1)**(-1) * Q(2)**(-1) * cos(Q(3))*sin(Q(3))**(-2)
Vep = Vep -0.00014000206382_real64 * Q(1)**(-2)
Vep = Vep -0.00007000103191_real64 * Q(1)**(-2) * cos(Q(3))**(2)*sin(Q(3))**(-2)
Vep = Vep -0.00007201429494_real64 * Q(2)**(-2)
Vep = Vep -0.00003600714747_real64 * Q(2)**(-2) * cos(Q(3))**(2)*sin(Q(3))**(-2)

END SUBROUTINE Tana_F2_F1_Vep
