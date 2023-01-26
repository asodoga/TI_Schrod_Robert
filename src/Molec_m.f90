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


!===============================================================================
module Molec_m
  USE NumParameters_m
  USE Tana_m

  IMPLICIT NONE
  !Real(kind=Rk) :: mass = HALF
  Real(kind=Rk) :: mass3(3) != HALF
  Real(kind=Rk) :: massr3(3)
  !Real(kind=Rk) :: mass = 1744.44536_RK
  Real(kind=Rk) :: mass = ONE
  Real(kind=Rk) :: De = 0.2250_RK
  Real(kind=Rk) :: Re = 1.73290_RK
  Real(kind=Rk) :: alpha1 = 1.1741_RK
  Logical       :: QML=.True.

  TYPE Molec_t
    Integer                        :: ndim
    Real (kind=Rk)                 :: V0 ! le fond du puits
    Character(len=:), allocatable  :: sym_type !type of symmetry
    Character(len=:), allocatable  :: coord_type ! name of coord
    Character(len=:), allocatable  :: Model_type
  END TYPE Molec_t

  public :: Calc_pot,Set_mass,mass,mass3,massr3,Molec_t,Read_Molec,Calc_potsub,Tana_F2_F1_Vep

  contains

   SUBROUTINE Set_mass(massr3)
    Real(kind=Rk),intent(out)           :: massr3(3)
    Real(kind=Rk)                       :: mass3(3)
    !Logical,          parameter       :: debug = .true.
    Logical,         parameter         :: debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Set_mass'
      flush(out_unitp)
    END IF

    mass3(1) =  0.00100782503223_Rk*((6.02214076_Rk*TEN**23)*(9.1093837015_Rk*TEN**(-31)))**(-1)
    mass3(2) =  0.00100782503223_Rk*((6.02214076_Rk*TEN**23)*(9.1093837015_Rk*TEN**(-31)))**(-1)
    mass3(3) =  0.034968852682_Rk*((6.02214076_Rk*TEN**23)*(9.1093837015_Rk*TEN**(-31)))**(-1)


    massr3(1) = mass3(1)*mass3(3)*(mass3(1)+mass3(3))**(-1)
    massr3(2) = mass3(1)*mass3(3)*(mass3(1)+mass3(3))**(-1)
    massr3(3) = mass3(3)


    IF (debug) THEN
      write(out_unitp,*) 'END Set_mass'
      flush(out_unitp)
    END IF

  END SUBROUTINE Set_mass

  SUBROUTINE Read_Molec(Molec,ni)
  USE UtilLib_m
    Integer,             intent(in)     :: ni
    TYPE(Molec_t),       intent(inout)  :: Molec
    !Logical,             parameter     :: debug = .true.
    Logical,             parameter      :: debug = .false.
    Real(kind=Rk)                       :: V0
    Integer                             :: err_io
    Character (len=Name_len)            :: coord_type,sym_type,Model_type
    Integer                             :: ndim,ndimQ,nsurf,option,inb
    Logical                             :: adiabatic
    Character (len=16)                  :: pot_name
    Real(kind=Rk), allocatable          :: Mat_V(:,:)
    Real(kind=Rk), allocatable          :: H(:,:,:,:),QML0(:),Q0(:)
    Real(kind=Rk), allocatable          :: G(:,:,:)
    Real(kind=Rk), allocatable          :: F2(:,:),F1(:),ScaleQ(:),m(:)
    Real(kind=Rk)                       :: FR,Vep

    NAMELIST /pot/ coord_type,sym_type,Model_type,V0,ndim

    Model_type    = 'local'
    coord_type    = 'simple'
    sym_type      = 'sym'
    V0            = Zero
    ndim          = 1

    IF (debug) THEN
      write(out_unitp,*)'BEGINNING Read_Molec'
      write(out_unitp,*)'V0=',V0,'coord_type=',coord_type,'sym_type=',sym_type
      write(out_unitp,*)'Model_type=',Model_type
      flush(out_unitp)
    END IF

    Read(ni,nml=pot,IOSTAT=err_io)
    SELECT CASE (Model_type)
    CASE('QML')
      ndimQ      = 0
      nsurf     = 0
      pot_name  = 'read_model'
      adiabatic = .FALSE.
      option    = 0

      Call sub_Init_Qmodel(ndimQ,nsurf,pot_name,adiabatic,option)

      Molec%ndim  =  ndimQ

      Allocate(Mat_V(nsurf,nsurf))
      Allocate(H(nsurf,nsurf,ndimQ,ndimQ))
      Allocate(G(nsurf,nsurf,ndimQ))
      Allocate(F2(ndimQ,ndimQ))
      Allocate(F1(ndimQ))
      Allocate(Q0(ndimQ))
      Allocate(m(ndimQ))
      Allocate(QML0(ndimQ))
      Allocate(ScaleQ(ndimQ))

      Call get_Qmodel_Q0(QML0,option)

      SELECT CASE (sym_type)
      CASE ('ANTISYM')
       Q0(3)=QML0(1)
       Q0(1)=QML0(2)
       Q0(2)=QML0(2)

       Call Tana_F2_F1_Vep(F2,F1,Vep,Q0)

       Call sub_Qmodel_VGH(Mat_V,G,H,QML0)

       m(1)      = -one/(2*F2(1,1))
       ScaleQ(3) = sqrt(m(1)*sqrt(H(1,1,1,1)/m(1)))
       m(2)      = -one/(2*F2(2,2))
       Fr        = H(1,1,2,2)
       ScaleQ(1) = sqrt(m(2)*sqrt(HALF*Fr/(m(2))))
       ScaleQ(2) = sqrt(m(2)*sqrt(HALF*Fr/(m(2))))

       DO inb = 1, 3
         Write(*,*)inb, 'ScaleQ=',ScaleQ(inb)
       END DO

      CASE ('SYM')
       Q0(1)=QML0(1)
       Q0(2)=QML0(2)
       Q0(3)=QML0(3)

       Call Tana_F2_F1_Vep(F2,F1,Vep,Q0)
       Call sub_Qmodel_VGH(Mat_V,G,H,QML0)

       DO inb =1,3
         m(inb)      = -one/(2*F2(inb,inb))
         ScaleQ(inb) =sqrt(m(inb)*sqrt(H(1,1,inb,inb)/m(inb)))
         Write(*,*)inb, 'ScaleQ(inb)=',ScaleQ(inb)
       END DO

      CASE DEFAULT
          STOP 'sym_type is bad'
      END SELECT

    CASE('LOCAL')
      Molec%ndim  =  ndim
    CASE DEFAULT
      STOP 'Model_type is bad'
    END SELECT

    Molec%V0         = V0
    Molec%coord_type = coord_type
    Molec%sym_type   = sym_type
    Molec%Model_type = Model_type

    IF (debug) THEN
      write(out_unitp,*)'V0=',Molec%V0,'coord_type=',Molec%coord_type,'sym_type=',Molec%sym_type
      write(out_unitp,*)'Model_type=',Molec%Model_type
      write(out_unitp,*)'END Read_Molec'
      flush(out_unitp)
    END IF

    Deallocate(Mat_V)
    Deallocate(H)
    Deallocate(G)
    Deallocate(F2)
    Deallocate(F1)
    Deallocate(Q0)
    Deallocate(m)
    Deallocate(QML0)
    Deallocate(ScaleQ)

 END SUBROUTINE Read_Molec

 SUBROUTINE Calc_potsub(Calc_pot,Q,Molec)
 USE UtilLib_m
    Real(kind=Rk), intent(in)     :: Q(:)
    TYPE(Molec_t), intent(in)     :: Molec
    Real(kind=Rk), intent(inout)  :: Calc_pot
    Real(kind=Rk), allocatable    :: QQML(:)
    Real(kind=Rk), allocatable    :: Mat_V(:,:)
    Real(kind=Rk)                 :: A,B,Tmin,Tmax
    !Logical,       parameter     :: debug = .true.
    Logical,       parameter      :: debug = .false.
    Logical,       parameter      :: Pot_calc = .false.
    integer                       :: iq1,iq2,I


    IF (debug) THEN
      write(out_unitp,*)'BEGINNING Calc_potsub'
      write(out_unitp,*)
      write(out_unitp,*)
      flush(out_unitp)
    END IF

    SELECT CASE (Molec%Model_type)
    CASE ('QML')
      ALLocate(Mat_V(1,1))
      ALLocate(QQML(3))
      SELECT CASE (Molec%sym_type)
      CASE ('ANTISYM')
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
      CASE ('SYM')
       QQML(1)= Q(1)
       QQML(2)= Q(2)
       QQML(3)= Q(3)
      CASE DEFAULT
        STOP 'sym_type is bad'
      END SELECT
      Call sub_Qmodel_V(Mat_V,QQML)
      Calc_pot = Mat_V(1,1)
    CASE ('Local')
      Calc_pot = HALF * dot_product( Q,Q)! 0.5*x^2
      !Calc_pot =  dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! x^2+x^4
      !Calc_pot =  -TEN * dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! -10x^2+x^4
      !Calc_pot = HALF * dot_product( Q,Q)
      !Calc_pot = De*dot_product((1-exp(-alpha1*(Q-Re))),(1-exp(-alpha1*(Q-Re))))
      !Calc_pot = De*(1-exp(-alpha1*(Q(1)-Re)))**2
    END SELECT
!!!!!!!Calcul du potentiel!!!!
  IF (Pot_calc) THEN
    Open(1,file='Pot_r1_r2v.dat',status='replace')
    !I=0
    DO iq2 = 0,150
    DO iq1 = 0,150

     QQML(1)=   1.9_Rk+iq1*0.01_Rk
      QQML(2)=   1.9_Rk+iq2*0.01_Rk !2.46_Rk !+ iq*0.2_Rk
      QQML(3)=   1.64_rk !+ iq*0.5_Rk
      Call sub_Qmodel_V(Mat_V,QQML)
      Write(1,*) QQML,(Mat_V+460.59052688079885_Rk)*219474.631443_RK
    END DO
    Write(1,*)
    END DO
    Close(1)

    Open(2,file='Pot_r1_Theta.dat',status='replace')
    !I=0
    DO iq2 = 0,100
    DO iq1 = 0,150

      QQML(1) =   1.9_Rk+iq1*0.01_Rk
      QQML(2) =   2.46_Rk
      QQML(3) =   One + iq2*0.02_Rk
      Call sub_Qmodel_V(Mat_V,QQML)
      Write(2,*) QQML,(Mat_V+460.59052688079885_Rk)*219474.631443_RK
    END DO
    Write(2,*)
    END DO
    Close(2)
  END IF

    Deallocate(Mat_V)
    Deallocate(QQML)
    IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) Calc_pot
      write(out_unitp,*) 'END Calc_potsub'
      flush(out_unitp)
    END IF

 END SUBROUTINE Calc_potsub

 FUNCTION Calc_pot(Q)
   Real(kind=Rk)               :: Calc_pot
   Real(kind=Rk), intent(in)   :: Q(:)
   Real(kind=Rk), allocatable  :: QQML(:)
   Real(kind=Rk), allocatable  :: Mat_V(:,:)
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

      Call sub_Qmodel_V(Mat_V,QQML)
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
   Deallocate(Mat_V)
   Deallocate(QQML)
  END FUNCTION Calc_pot

End module Molec_m
