!===============================================================================
!  This file is part of TI_Schrod_Robert.
!
!  TI_Schrod_Robert is a free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!  TI_Schrod_Robert is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with TI_Schrod_Robert.  If not, see <http://www.gnu.org/licenses/>.
!
!   Copyright 2022 Robert AFANSOUNOUDJI [1]

!     with contributions of:
!     Komi SODOGA      [1]
!     David Lauvergnat [2]
![1]:Laboratoire de Physique des Matériaux et des Composants
!à Semi-Conducteurs (LPMCS), UNIVERSITÉ DE LOME
![2]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France

!===============================================================================
module Molec_m
  USE NumParameters_m
  implicit none
  private

  !real(kind=Rk) :: mass = HALF
  real(kind=Rk) :: mass3(3) != HALF
  real(kind=Rk) :: massr3(3)

  !real(kind=Rk) :: mass = 1744.44536_RK
  real(kind=Rk) :: mass = ONE
  real(kind=Rk) :: De = 0.2250_RK
  real(kind=Rk) :: Re = 1.73290_RK
  real(kind=Rk) :: alpha1 = 1.1741_RK
  logical       :: QML=.True.

  public :: Calc_pot,Set_mass,mass,mass3


  contains

   SUBROUTINE Set_mass(massr3)

    real(kind=Rk),intent(out)           :: massr3(3)
    real(kind=Rk)                       :: mass3(3)
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.


    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Set_mass'
      flush(out_unitp)
    END IF

    mass3(1) =  0.00100782503223_Rk*((6.02214076_Rk*TEN**23)*(9.1093837015_Rk*TEN**(-31)))**(-1)
    mass3(2) =  0.00100782503223_Rk*((6.02214076_Rk*TEN**23)*(9.1093837015_Rk*TEN**(-31)))**(-1)
    mass3(3) =  0.034968852682_Rk*((6.02214076_Rk*TEN**23)*(9.1093837015_Rk*TEN**(-31)))**(-1)


    !write(6,*) 'mass',mass3(:)


    massr3(1) = mass3(1)*mass3(3)*(mass3(1)+mass3(3))**(-1)
    massr3(2) = mass3(1)*mass3(3)*(mass3(1)+mass3(3))**(-1)
    massr3(3) =  mass3(3)


    IF (debug) THEN
      write(out_unitp,*) 'END Set_mass'
      flush(out_unitp)
    END IF

  END SUBROUTINE Set_mass



  FUNCTION Calc_pot(Q)
   real(kind=Rk)               :: Calc_pot
   real(kind=Rk), intent(in)   :: Q(:)
   REAL(kind=Rk), ALLOCATABLE  :: QQML(:)
   REAL(kind=Rk), ALLOCATABLE  :: Mat_V(:,:)
    !  write(*,*) "Q",Q
   IF(QML) THEN
     ALLocate(Mat_V(1,1))
     ALLocate(QQML(3))
    ! QQML(1) = a          (Radian)
    ! QQML(2) = 1/2(R1+R2) (Bohr)
    ! QQML(3) = 1/2(R1-R2) (Bohr)
    !Q(1)==R1
    !Q(2)==R2
    !Q(3)==a
    ! QQML(:)= [Q(3),HALF*(Q(1)+Q(2)),HALF*(Q(1)-Q(2))]

      QQML(1)= Q(3)
      QQML(2)= HALF*(Q(1)+Q(2))
      QQML(3)= HALF*(Q(1)-Q(2))

      !write(*,*) "QQML",QQML

      CALL sub_Qmodel_V(Mat_V,QQML)
      !write(*,*) "Mat_V",Mat_V
      Calc_pot = Mat_V(1,1)
   ELSE
      Calc_pot = HALF * dot_product( Q,Q)! 0.5*x^2
    !Calc_pot =  dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! x^2+x^4
    !Calc_pot =  -TEN * dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! -10x^2+x^4
    !Calc_pot = HALF * dot_product( Q,Q)
    !Calc_pot = De*dot_product((1-exp(-alpha1*(Q-Re))),(1-exp(-alpha1*(Q-Re))))
   ! Calc_pot = De*(1-exp(-alpha1*(Q(1)-Re)))**2
   END IF
  END FUNCTION Calc_pot

end module Molec_m
