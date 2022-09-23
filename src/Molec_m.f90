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
  USE Tana_m

  IMPLICIT NONE


  !real(kind=Rk) :: mass = HALF
  real(kind=Rk) :: mass3(3) != HALF
  real(kind=Rk) :: massr3(3)

  !real(kind=Rk) :: mass = 1744.44536_RK
  real(kind=Rk) :: mass = ONE
  real(kind=Rk) :: De = 0.2250_RK
  real(kind=Rk) :: Re = 1.73290_RK
  real(kind=Rk) :: alpha1 = 1.1741_RK
  logical       :: QML=.True.

  TYPE Molec_t
    integer                        :: ndim
    real (kind=Rk)                 :: V0 ! le fond du puits
    character(len=:), allocatable  :: sym_type !type of symmetry
    character(len=:), allocatable  :: coord_type ! name of coord
    character(len=:), allocatable  :: Model_type
  END TYPE Molec_t

  public :: Calc_pot,Set_mass,mass,mass3,massr3,Molec_t,Read_Molec,Calc_potsub,Tana_F2_F1_Vep
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
    massr3(3) = mass3(3)


    IF (debug) THEN
      write(out_unitp,*) 'END Set_mass'
      flush(out_unitp)
    END IF

  END SUBROUTINE Set_mass

  SUBROUTINE Read_Molec(Molec,ni)
  USE UtilLib_m
    integer,             intent(in)     :: ni
    TYPE(Molec_t),       intent(inout)  :: Molec
    !logical,             parameter      :: debug = .true.
    logical,             parameter      :: debug = .false.
    REAL(kind=Rk)                       :: V0
    integer                             :: err_io
    character (len=Name_len)            :: coord_type,sym_type,Model_type
    integer                             :: ndim,ndimQ,nsurf,option,inb
    logical                             :: adiabatic
    character (len=16)                  :: pot_name
    REAL(kind=Rk), ALLOCATABLE          :: Mat_V(:,:)
    REAL(kind=Rk), ALLOCATABLE          :: H(:,:,:,:),QML0(:),Q0(:)
    REAL(kind=Rk), ALLOCATABLE          :: G(:,:,:)
    REAL(kind=Rk), ALLOCATABLE          :: F2(:,:),F1(:),ScaleQ(:),m(:)
    REAL(kind=Rk)                       :: FR,Vep

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
      CALL sub_Init_Qmodel(ndimQ,nsurf,pot_name,adiabatic,option)

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

      CALL get_Qmodel_Q0(QML0,option)

      SELECT CASE (sym_type)
      CASE ('ANTISYM')
       Q0(3)=QML0(1)
       Q0(1)=QML0(2)
       Q0(2)=QML0(2)
       Call Tana_F2_F1_Vep(F2,F1,Vep,Q0)
       CALL sub_Qmodel_VGH(Mat_V,G,H,QML0)
       m(1) = -one/(2*F2(1,1))
       ScaleQ(3)=sqrt(m(1)*sqrt(H(1,1,1,1)/m(1)))
       m(2) = -one/(2*F2(2,2))
       Fr=H(1,1,2,2)
       ScaleQ(1)= sqrt(m(2)*sqrt(HALF*Fr/(m(2))))
       ScaleQ(2)= sqrt(m(2)*sqrt(HALF*Fr/(m(2))))
       DO inb =1,3
         Write(*,*)inb, 'ScaleQ=',ScaleQ(inb)
       END DO
      CASE ('SYM')
       Q0(1)=QML0(1)
       Q0(2)=QML0(2)
       Q0(3)=QML0(3)
       Call Tana_F2_F1_Vep(F2,F1,Vep,Q0)
       CALL sub_Qmodel_VGH(Mat_V,G,H,QML0)
       DO inb =1,3
         m(inb) = -one/(2*F2(inb,inb))
         ScaleQ(inb)=sqrt(m(inb)*sqrt(H(1,1,inb,inb)/m(inb)))
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

 END SUBROUTINE Read_Molec

 SUBROUTINE Calc_potsub(Calc_pot,Q,Molec)
  USE UtilLib_m
    REAL(kind=Rk),       intent(in)     :: Q(:)
    TYPE(Molec_t),       intent(in)     :: Molec
    REAL(kind=Rk),       intent(inout)  :: Calc_pot
    REAL(kind=Rk), ALLOCATABLE          :: QQML(:)
    REAL(kind=Rk), ALLOCATABLE          :: Mat_V(:,:)
    !logical,             parameter      :: debug = .true.
    logical,             parameter      ::debug = .false.


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
      CASE ('sym')
       QQML(1)= Q(1)
       QQML(2)= Q(2)
       QQML(3)= Q(3)
      CASE DEFAULT
          STOP 'sym_type is bad'
      END SELECT
      !write(*,*) "QQML",QQML
      CALL sub_Qmodel_V(Mat_V,QQML)
      !write(*,*) "Mat_V",Mat_V
      Calc_pot = Mat_V(1,1)
    CASE ('local')
      Calc_pot = HALF * dot_product( Q,Q)! 0.5*x^2
     !Calc_pot =  dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! x^2+x^4
     !Calc_pot =  -TEN * dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! -10x^2+x^4
     !Calc_pot = HALF * dot_product( Q,Q)
     !Calc_pot = De*dot_product((1-exp(-alpha1*(Q-Re))),(1-exp(-alpha1*(Q-Re))))
     ! Calc_pot = De*(1-exp(-alpha1*(Q(1)-Re)))**2
    END SELECT
    IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) Calc_pot
      write(out_unitp,*) ' END Calc_potsub'
      flush(out_unitp)
    END IF
 END SUBROUTINE Calc_potsub

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
