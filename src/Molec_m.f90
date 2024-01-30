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
  USE Pot_m

  IMPLICIT NONE
  !TYPE(Para_t)     :: Para
  !TYPE(Cage_H2O_t) :: Cage_H2O
  !Real(kind=Rk)   :: mass = HALF
  Real(kind=Rk)    :: mass3(3) != HALF
  Real(kind=Rk)    :: massr3(3)
  !Real(kind=Rk)   :: mass = 1744.44536_RK
  Real(kind=Rk)    :: mass = ONE
  Real(kind=Rk)    :: De = 0.2250_RK
  Real(kind=Rk)    :: Re = 1.73290_RK
  Real(kind=Rk)    :: alpha1 = 1.1741_RK
  Logical          :: QML=.True.


  TYPE Molec_t
    Integer                        :: ndim
    Real (kind=Rk)                 :: V0 ! le fond du puits
    Character(len=:), allocatable  :: sym_type !type of symmetry
    Character(len=:), allocatable  :: coord_type ! name of coord
    Character(len=:), allocatable  :: Model_type
    TYPE(Para_CO_nWater_t)         :: Para_CO_nWater
    TYPE(Cage_mole_t)              :: Cage_mole
  END TYPE Molec_t

  public :: mass,mass3,massr3
  public :: Molec_t,Read_Molec,Calc_potsub,Tana_F2_F1_Vep

  contains


  SUBROUTINE Read_Molec(Molec,ni)
  USE UtilLib_m
    Integer,             intent(in)     :: ni
    TYPE(Molec_t),       intent(inout)  :: Molec
    !TYPE(Para_t)                        :: Para
    !TYPE(Cage_H2O_t)                    :: Cage_H2O
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

    Model_type    = 'LOCAL'
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

      ndimQ     = 0 !H2cl+
      !ndimQ     = 5
    !  nsurf     = 1
      nsurf     = 0 ! H2cl+
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
       Q0(3) = QML0(1)
       Q0(1) = QML0(2)
       Q0(2) = QML0(2)

       Call Tana_F2_F1_Vep(F2,F1,Vep,Q0)

       Call sub_Qmodel_VGH(Mat_V,G,H,QML0)

       m(1)      = -one/(2*F2(1,1))
       ScaleQ(3) = sqrt(m(1)*sqrt(H(1,1,1,1)/m(1)))
       m(2)      = -one/(2*F2(2,2))
       Fr        = H(1,1,2,2)
       ScaleQ(1) = sqrt(m(2)*sqrt(HALF*Fr/(m(2))))
       ScaleQ(2) = sqrt(m(2)*sqrt(HALF*Fr/(m(2))))

       DO inb = 1, Molec%ndim
         Write(*,*)inb, 'ScaleQ=',ScaleQ(inb)
       END DO

      CASE ('SYM')

        DO inb = 1,Molec%ndim
          Q0(inb) = QML0(inb)
        END DO


       Call Tana_F2_F1_Vep(F2,F1,Vep,Q0)
       Call sub_Qmodel_VGH(Mat_V,G,H,QML0)

       DO inb = 1,Molec%ndim
         m(inb)      = -one/(2*F2(inb,inb))
         ScaleQ(inb) =  sqrt(m(inb)*sqrt(H(1,1,inb,inb)/m(inb)))
         Write(*,*)inb, 'ScaleQ(inb)=',ScaleQ(inb)
       END DO
       Deallocate(Mat_V)
       Deallocate(H)
       Deallocate(G)
       Deallocate(F2)
       Deallocate(F1)
       Deallocate(Q0)
       Deallocate(m)
       Deallocate(QML0)
       Deallocate(ScaleQ)
      CASE DEFAULT
          STOP 'sym_type is bad'
      END SELECT
    CASE('RPotlib')
      Call Set_para(Molec%Para_CO_nWater)
      Call Read_Cage(Molec%Cage_mole,Molec%Para_CO_nWater)

    CASE('LOCAL')
      Molec%ndim  =  ndim
    CASE DEFAULT
      STOP 'Model_type is bad in Read_Molec'
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
    Real(kind=Rk), intent(in)     :: Q(:)
    TYPE(Molec_t), intent(in)     :: Molec
    !TYPE(Para_t)                  :: Para
    !TYPE(Cage_H2O_t)              :: Cage_H2O
    Real(kind=Rk), intent(inout)  :: Calc_pot
    Real(kind=Rk), allocatable    :: QQML(:)
    Real(kind=Rk)                 :: Q2(6)
    Real(kind=Rk), allocatable    :: Mat_V(:,:)
    Real(kind=Rk)                 :: A,B,Tmin,Tmax
    !Logical,       parameter     :: debug = .true.
    Logical,       parameter      :: debug = .false.
    Logical,       parameter      :: Pot_calc = .false.
    integer                       :: iq1,iq2,I, inb


    IF (debug) THEN
      write(out_unitp,*)'BEGINNING Calc_potsub'
      write(out_unitp,*)
      write(out_unitp,*)
      flush(out_unitp)
    END IF

    SELECT CASE (Molec%Model_type)
    CASE ('QML')
      ALLocate(Mat_V(1,1))
      ALLocate(QQML(Molec%ndim))
      SELECT CASE (Molec%sym_type)
      CASE ('ANTISYM')
        QQML(:)= [Q(3),HALF*(Q(1)+Q(2)),HALF*(Q(1)-Q(2))]
      CASE ('SYM')

       DO inb = 1,Molec%ndim
         QQML(inb)= Q(inb)
       END DO

      CASE DEFAULT
        STOP 'sym_type is bad '
      END SELECT
      Call sub_Qmodel_V(Mat_V,QQML)
      Calc_pot = Mat_V(1,1)
      Deallocate(Mat_V)
      Deallocate(QQML)
      write(*,*) 'Calc_pot=',Calc_pot
    CASE('RPotlib')
    !  write(*,*) 'size(Q)=',size(Q2)
      !Q2(:) = one
      Call Pot_CO_Cage( Calc_pot,Molec%Para_CO_nWater,Molec%Cage_mole,Q)
      !Write(*,*) 'Calc_pot=',Calc_pot
    !  write(*,*) 'size(Q)=',size(Q2)
      !Stop 'Robert'
    CASE ('LOCAL')
    !  Stop 'Robert'
      Calc_pot = HALF * dot_product( Q,Q)! 0.5*x^2
      !write(*,*) 'Calc_pot=',Calc_pot
    CASE DEFAULT
      STOP 'Model_type is bad in Calc_potsub'
    END SELECT
!!!!!!!Calcul du potentiel!!!!
  IF (Pot_calc) THEN
    ALLocate(Mat_V(1,1))
    ALLocate(QQML(3))
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
    Deallocate(Mat_V)
    Deallocate(QQML)
  END IF


  IF (debug) THEN
    write(out_unitp,*)
    write(out_unitp,*) Calc_pot
    write(out_unitp,*) 'END Calc_potsub'
    flush(out_unitp)
  END IF

 END SUBROUTINE Calc_potsub


END MODULE  Molec_m
