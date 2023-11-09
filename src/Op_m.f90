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
!This module allows to do the action of the Hamiltonian operator on the wave
! packet and then builds the Hamiltonian matrixAllocate(Basis%tab_Smolyak(I)%NDindexb%Li_smol(Basis%NDindexl%NDim))
!====================================================
MODULE Op_m
  USE NumParameters_m
  USE diago_m
  USE NDindex_m
  USE Molec_m
  USE Basis_m, only : Basis_t
  USE psi_m

  IMPLICIT  NONE

  TYPE :: Grid_t
    Real(kind=Rk),    allocatable :: Vec(:)
    Integer                       :: DerivIndex(2)
  END TYPE Grid_t

  TYPE :: Op_t
    Real(kind=Rk),   allocatable  :: Vg(:)
    Real(kind=Rk),   allocatable  :: Vb(:)
    TYPE (Grid_t),   allocatable  :: Grid(:)
    TYPE (Basis_t),   pointer     :: Basis => Null()
    TYPE (Molec_t),   pointer     :: Molec => Null()
    Real (kind=Rk),   allocatable :: Scalar_g(:) ! The scalar part of the operator (The  for Hamiltonian)
    Real (kind=Rk),   allocatable :: RMat(:,:) ! The matrix of the Hamiltonian operator
    TYPE (Op_t),      allocatable :: Op_smol(:)
  END TYPE Op_t

  Public :: Op_t,Write_Op,Set_Op,dealloc_Op,calc_OpPsi,Diago_Op
  Public :: OpPsi_gridnD,Grid_t,Make_mat,Local_F2_F1_old

CONTAINS

  SUBROUTINE alloc_Op(Op,nb)
    TYPE(Op_t),  intent(inout) :: Op
    Integer,     intent(in)    :: nb
    Call dealloc_Op(Op)
    IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'
    Allocate(Op%RMat(nb,nb))
  END SUBROUTINE alloc_Op

  SUBROUTINE dealloc_Op(Op)
    TYPE(Op_t), intent(inout) :: Op
    nullify(Op%Basis)
    IF (Allocated(Op%RMat)) THEN
      Deallocate(Op%RMat)
    END IF
    IF (Allocated(Op%Scalar_g)) THEN
      Deallocate(Op%Scalar_g)
    END IF
  END SUBROUTINE dealloc_Op

  SUBROUTINE Write_Op(Op)
   USE UtilLib_m, only : Write_RMat
    TYPE(Op_t), intent(in) :: Op
    IF (associated(Op%Basis)) THEN
      Write(out_unitp,*) ' The basis is linked to Op.'
    ELSE
      Write(out_unitp,*) ' The basis is NOT linked to Op.'
    END IF

    IF (Allocated(Op%RMat)) THEN
      Write(out_unitp,*) 'Writing Op (Real):'
      Write(out_unitp,*)
      Call Write_RMat(Op%RMat,out_unitp,5,name_info='Op%Rmat')
      Write(out_unitp,*) 'END Writing Op'
    END IF

  END SUBROUTINE Write_Op

  SUBROUTINE OpPsi_basis_DP(OpPsi_B,Psi_B,Op)
   USE Basis_m
   USE Molec_m
    IMPLICIT  NONE
    TYPE (Op_t),     intent(inout)      :: Op
    Real (kind=Rk), intent(out)         :: OpPsi_B(:)
    Real (kind=Rk), intent(in)          :: Psi_B(:)
    !Logical,          parameter        :: debug = .true.
    Logical,         parameter          :: debug = .false.
    Real (kind=Rk), allocatable         :: Psi_G(:)
    Real (kind=Rk), allocatable         :: OpPsi_G(:)

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING OpPsi_basis_DP'
     Call Write_basis(Op%Basis)
     flush(out_unitp)
    END IF
    Allocate(Psi_G(Op%Basis%nq))
    Allocate(OpPsi_G(Op%Basis%nq))

    CALL BasisTOGrid_Basis_rapide1(Psi_G,Psi_B,Op%Basis)

    CALL OpPsi_grid(OpPsi_G,Psi_G,Op)
    CALL GridTOBasis_Basis_rapide1( OpPsi_B,OpPsi_G,Op%Basis)
    Deallocate(Psi_G)
    Deallocate(OpPsi_G)
    IF (debug) THEN
     Write(out_unitp,*) 'END OpPsi_basis_DP'
     flush(out_unitp)
    END IF
  END SUBROUTINE OpPsi_basis_DP

  SUBROUTINE OpPsi_basis_Smol(OpPsi_Bc,Psi_Bc,Op)
   USE Basis_m
   USE Molec_m
    IMPLICIT  NONE
    TYPE (Op_t),     intent(inout)      :: Op
    Real (kind=Rk),  intent(in)         :: Psi_Bc(:)
    Real (kind=Rk),  intent(out)        :: OpPsi_Bc(:)
    Integer                             :: I_smol,I_compt,n, I_smol1
    Real (kind=Rk), allocatable         :: Psi_B(:)
    Real (kind=Rk), allocatable         :: OpPsi_B(:)

    OpPsi_Bc = ZERO

    DO I_smol = 1,Op%Basis%NDindexl%Nterm
     Allocate(Psi_B(Op%Basis%tab_Smolyak(I_smol)%nb))
     Allocate(OpPsi_B(Op%Basis%tab_Smolyak(I_smol)%nb))

     DO n = 1,Op%Basis%tab_Smolyak(I_smol)%nb
      I_compt = Op%Basis%Tab_Smolyak(I_smol)%Ind_map(n)
      Psi_B(n) = Psi_Bc(I_compt)
     END DO

     CALL OpPsi_basis_DP(OpPsi_B,Psi_B,Op%Op_Smol(I_smol))

     DO n = 1,Op%Basis%tab_Smolyak(I_smol)%nb
        I_compt = Op%Basis%Tab_Smolyak(I_smol)%Ind_map(n)
        OpPsi_Bc(I_compt) = OpPsi_Bc(I_compt) + Op%Basis%Smolyak%DSmol(I_smol)*OpPsi_B(n)
     END DO


     Deallocate(Psi_B)
     Deallocate(OpPsi_B)
    END DO

  END SUBROUTINE OpPsi_basis_Smol

  SUBROUTINE OpPsi_basis(OpPsi_B,Psi_B,Op)
   USE Basis_m
   USE Molec_m
    IMPLICIT  NONE
    TYPE (Op_t),     intent(inout)  :: Op
    Real (kind=Rk),  intent(in)     :: Psi_B(:)
    Real (kind=Rk),  intent(out)    :: OpPsi_B(:)

    SELECT CASE (Op%Basis%Basis_name)
     CASE ('dp')
      CALL OpPsi_basis_DP(OpPsi_B,Psi_B,Op)
     CASE ('smolyak')
      CALL OpPsi_basis_Smol(OpPsi_B,Psi_B,Op)
     CASE default
      STOP 'ERROR in Read_Basis: no default basis.'
    END SELECT

  END SUBROUTINE OpPsi_basis

  SUBROUTINE Make_mat(Op)
   USE Basis_m
   USE Molec_m
    IMPLICIT  NONE
    TYPE (Op_t),  intent(inout)  :: Op
    Integer                      :: ib
    Real (kind=Rk), allocatable  :: Psi_B(:)

    Allocate(Op%RMat(Op%Basis%nb,Op%Basis%nb))
    Allocate(Psi_B(Op%Basis%nb))
    DO ib=1,Op%Basis%nb
      Psi_B(:) = ZERO ; Psi_B(ib) = ONE
      CALL OpPsi_basis(Op%RMat(:,ib),Psi_B,Op)
    END DO
    Deallocate(Psi_B)

  END SUBROUTINE Make_mat

  SUBROUTINE  Local_F2_F1_old(F2,F1,R,m1,m2,Ndim,Q)
  USE NumParameters_m
  USE UtilLib_m
   IMPLICIT NONE
   Integer, intent(in)           :: Ndim
   Real (kind=Rk), intent(in)    :: Q(:)
   Real (kind=Rk), intent(inout) :: F1(:)
   Real (kind=Rk), intent(inout) :: F2(:,:)
   Real (kind=Rk), intent(in)    :: R,m1,m2
   Real (kind=Rk)                :: mr,mt
   !Logical,          parameter   :: debug = .true.
   Logical,         parameter   :: debug = .false.
   Integer                       :: i,j

   IF (debug) THEN
    Write(out_unitp,*) 'BEGINNING Local_F2_F1_old'
    flush(out_unitp)
   END IF

    !Write(*,*) Q
    mt = m1+m2
    mr = (m1*m2)/mt

    F2 = Zero
    F1 = Zero

    DO i=1,Ndim-2
      F2(i,i) = - HALF/mt ! F2 for x,y,z
    END DO

    F2(Ndim-1,Ndim-1) = -HALF/(mr*R**2) !F2 for theta
!Write(*,*) (sin(Q(Ndim-1)))**2,Q(Ndim-1)
    F2(Ndim,Ndim)     = -HALF/((mr*R**2)*(sin(Q(Ndim-1)))**2) ! f2 for phi
!STOP 'RObert Aa'
    F1(Ndim-1)        = -HALF/((mr*R*R*tan(Q(Ndim-1)))) !F1 for theta

    IF (debug) THEN
      Do i=1,Ndim
       Write(out_unitp,*)'F2',F2(i,i)
        Write(out_unitp,*)'F1',F1(i)
      END DO

     Write(out_unitp,*) 'END Local_F2_F1_old'
     flush(out_unitp)
    END IF
  END SUBROUTINE Local_F2_F1_old

  SUBROUTINE  Local_F2_F1_6D(F2,F1,R,m1,m2,Ndim,Q)
  USE NumParameters_m
  USE UtilLib_m
   IMPLICIT NONE
   Integer, intent(in)           :: Ndim
   Real (kind=Rk), intent(in)    :: Q(:)
   Real (kind=Rk), intent(inout) :: F1(:)
   Real (kind=Rk), intent(inout) :: F2(:,:)
   Real (kind=Rk), intent(in)    :: R,m1,m2
   Real (kind=Rk)                :: mr,mt
   Logical,          parameter   :: debug = .true.
   !Logical,         parameter   :: debug = .false.
   Integer                       :: i,j

   IF (debug) THEN
    Write(out_unitp,*) 'BEGINNING Local_F2_F1_6D'
    flush(out_unitp)
   END IF

    mt = m1+m2
    mr = (m1*m2)/mt

    F2 = Zero
    F1 = Zero

    DO i=1,Ndim-3
      F2(i,i) = - HALF/mt ! F2 for x,y,z
    END DO

    F2(Ndim-1,Ndim-1) = - HALF/(mr*Q(Ndim-2)**2) !F2 for theta

    F2(Ndim,Ndim)     = - HALF/((mr*Q(Ndim-2)**2)*(sin(Q(Ndim-1)))**(2)) ! f2 for phi

    F1(Ndim-1)        = - HALF/((mr*Q(Ndim-2)**2*tan(Q(Ndim-1)))) !F1 for theta

    IF (debug) THEN
      Do i=1,Ndim
       Write(out_unitp,*)'F2',F2(i,i)
        Write(out_unitp,*)'F1',F1(i)
      END DO

     Write(out_unitp,*) 'END Local_F2_F1_6D'
     flush(out_unitp)
    END IF
  END SUBROUTINE Local_F2_F1_6D


  SUBROUTINE Set_Op(Op,Basis,Molec)
   USE Basis_m
   USE Molec_m
   USE NDindex_m
    IMPLICIT  NONE
    TYPE(Op_t),     intent(inout)       :: Op
    TYPE (Basis_t), intent(in),  target :: Basis
    TYPE (Molec_t), intent(in),  target :: Molec
    !Logical,          parameter         :: debug = .true.
    Logical,         parameter          :: debug = .false.
    Integer                             :: ib,iq
    Logical                             :: Simple = .false.

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING Set_Op'
     Call Write_basis(Basis)
     flush(out_unitp)
    END IF

    IF (Simple) THEN
     Op%Basis => Basis  !Initialization of Op%Basis by the values ​​of Basis
     Op%Molec => Molec

     Allocate(Op%Grid(3))
     Allocate(Op%Grid(1)%Vec(Op%Basis%nq))
     Allocate(Op%Grid(2)%Vec(Op%Basis%nq))
     Allocate(Op%Grid(3)%Vec(Op%Basis%nq))

     Op%Grid(1)%DerivIndex(1) = 1
     Op%Grid(1)%DerivIndex(2) = 1
     Op%Grid(2)%DerivIndex(1) = 2
     Op%Grid(2)%DerivIndex(2) = 2
     Op%Grid(3)%DerivIndex(1) = 3
     Op%Grid(3)%DerivIndex(2) = 3

     Op%Grid(1)%Vec(:) = -HALF/mass
     Op%Grid(2)%Vec(:) = -HALF/mass
     Op%Grid(3)%Vec(:) = -HALF/mass
    ELSE
     Call Set_grid_Op(Op,Basis,Molec)
    END IF

    IF (debug) THEN
     Write(out_unitp,*) 'END Set_Op'
     flush(out_unitp)
    END IF
  END SUBROUTINE Set_Op

  SUBROUTINE Set_grid_Op(Op,Basis,Molec)
   USE Basis_m
   USE NDindex_m

    IMPLICIT  NONE
    TYPE(Op_t),     intent(inout)       :: Op
    TYPE (Basis_t), intent(in),  target :: Basis
    TYPE (Molec_t), intent(in),  target :: Molec
    Real(kind=Rk)                       :: F,Vep,V
    Real(kind=Rk), allocatable          :: Q(:),F2(:,:),F1(:)
    Real(kind=Rk)                       :: massCO(2)
    !Logical,          parameter         :: debug = .true.
    Logical, parameter                  :: debug = .false.
    Logical                             :: Endloop_q
    Logical                             :: Endloop_b
    Integer, allocatable                :: tab_iq(:)
    Integer                             :: ib,iq,i1,i2,inb,iterm,Nterm,ndim
    Integer                             :: inq1,inq2,I_smol

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING Set_grid_Op'
     !Call Write_basis(Op%Basis)
     flush(out_unitp)
    END IF

    SELECT CASE (Basis%Basis_name)
     CASE ('dp')
      Op%Basis => Basis  !Initialization of Op%Basis by the values ​​of Basis
      Op%Molec => Molec

      Allocate(Q(size(Op%Basis%tab_basis)))
      Allocate(F1(size(Op%Basis%tab_basis)))
      Allocate(F2(size(Op%Basis%tab_basis),size(Op%Basis%tab_basis)))
      Nterm  =  (Op%Basis%NDindexq%Ndim + 2)*(Op%Basis%NDindexq%Ndim + 1)/2
      Allocate(Op%Grid(Nterm))

      DO iterm = 1,size(Op%Grid)
       Allocate(Op%Grid(iterm)%Vec(Op%Basis%nq))
      END DO
      iterm = 1
      Op%Grid(iterm)%DerivIndex(:)  = [0,0]
      DO i2 = 1,size(Op%Basis%tab_basis)
      DO i1 = i2,size(Op%Basis%tab_basis)
       iterm = iterm + 1
       Op%Grid(iterm)%DerivIndex(:)  = [i2,i1]
      END DO
      END DO

      DO i1=1,size(Op%Basis%tab_basis)
       iterm = iterm+1
       Op%Grid(iterm)%DerivIndex(:)  = [i1,0]
      ! Write(out_unitp,*) iterm, Op%Grid(iterm)%DerivIndex(:)
      END DO
      Allocate(Tab_iq(size(Op%Basis%tab_basis)))

      Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)

      Iq=0
      DO
       Iq=Iq+1

       Call increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)
       IF (Endloop_q) exit

       DO inb = 1,size(Op%Basis%tab_basis)
         Q(inb) =  Op%Basis%tab_basis(inb)%x(tab_iq(inb),1)
       END DO
       SELECT CASE (Op%Molec%Model_type)
       CASE('QML')
         Call Tana_F2_F1_Vep(F2,F1,Vep,Q)
       CASE('RPotlib')

         massCO(1)   =  12.0000000_Rk
         massCO(2)   =  15.99491461957_Rk
         Vep=Zero
         Call Local_F2_F1_old(F2,F1,Op%Molec%Para_CO_nWater%ReCO,&
         massCO(1),massCO(2),Op%Basis%NDindexq%Ndim,Q)

       CASE DEFAULT
         STOP 'Model_type is bad'
       END SELECT

       Call Calc_potsub(V,Q,Op%Molec)
!Write(*,*)size(Op%Grid)
!stop 'CCC'
       DO iterm = 1,size(Op%Grid)
        inq1  = Op%Grid(iterm)%DerivIndex(1)
        inq2  = Op%Grid(iterm)%DerivIndex(2)
        IF (inq1 > 0 .and. inq2 > 0 )THEN
          Op%Grid(iterm)%Vec(Iq)  =  F2(inq1,inq2)
        ELSE IF (inq1 > 0 .and. inq2 == 0) THEN
         Op%Grid(iterm)%Vec(Iq)  =  F1(inq1)
        ELSE IF (inq1 == 0 .and. inq2 == 0) THEN
         Op%Grid(iterm)%Vec(Iq)     =   V+Vep
        END IF
       END DO
      END DO

      Deallocate(Tab_iq)
     CASE ('smolyak')

      Op%Molec => Molec
      Op%Basis => Basis
      Allocate(Op%Op_smol(Basis%NDindexl%Nterm))

      DO I_smol=1, Basis%NDindexl%Nterm
        Op%Op_Smol(I_smol)%Basis => Basis%Tab_Smolyak(I_smol)
        CALL Set_Grid_Op(Op%Op_Smol(I_smol),Basis%Tab_Smolyak(I_smol),Op%Molec)
      END DO

     CASE default
       STOP 'SSERROR in Read_Basis: no default basis.'
    END SELECT

    IF (debug) THEN
      Write(out_unitp,*) 'END Set_drid_Op'
      flush(out_unitp)
    END IF
  END SUBROUTINE Set_grid_op


  SUBROUTINE Diago_Op(Op)
   USE Molec_m
   USE Basis_m
    TYPE(Op_t),     intent(inout)       :: Op
    !Logical,          parameter         :: debug = .true.
    Logical,         parameter          :: debug = .false.
    Integer                             :: jb,inb,ndim
    Integer                             :: ib,iq,L
    Logical                             :: Endloop_q
    Integer,        allocatable         :: tab_iq(:),NDinit(:),NDend(:)
    Real (kind=Rk), allocatable         :: EigenVal(:), EigenVec(:,:), Psi_g(:), DifEigen(:)
    Real (kind=Rk), allocatable         :: Qmoy(:), Rho_r(:,:), omega(:)
    Real (kind=Rk)                      :: W,Wr

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING Diago_Op'
     flush(out_unitp)
    END IF

    Allocate(DifEigen(Op%Basis%nb))
    !Allocate(omega(Op%Basis%nb))
    Allocate(EigenVal(Op%Basis%nb))
    Allocate(EigenVec(Op%Basis%nb,Op%Basis%nb))
    !Call Diago_Arpack(psi,Ene,nb_diago,max_diago,max_it,para_H)
    !Call Diago_Arpack(EigenVec,EigenVal,Op%Basis%nb,Op%Basis%nb,Op%Basis%nb,Op%RMat)
    Call  diagonalization(Op%RMat,EigenVal,EigenVec,Op%Basis%nb)
    Write(out_unitp,*)
    Write(out_unitp,*)
    Write(out_unitp,*) 'eigenvalues = '
    Write(out_unitp,*) 'n','EigenVal','E_n-E_1','omega'
    DO ib=1,21!Op%Basis%nb
     DifEigen(ib)=(EigenVal(ib)-Op%Molec%V0)*219474.631443_RK-(EigenVal(1)-Op%Molec%V0)*219474.631443_RK
     Write(out_unitp,*) ib,(EigenVal(ib)-Op%Molec%V0)*219474.631443_RK,'DifEigen(ib)',DifEigen(ib)
    END DO

    Write(out_unitp,*)
    Write(out_unitp,*)

    IF (debug) THEN
     Allocate(Psi_g(Op%Basis%nq))
     Call BasisTOGrid_Basis(Psi_g, EigenVec(:,1),Op%Basis)

     Allocate(Rho_r(Op%Basis%nq,size(Op%Basis%tab_basis)))
     Allocate(Qmoy(size(Op%Basis%tab_basis)))
     Allocate(Tab_iq(size(Op%Basis%tab_basis)))
     Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
     Qmoy(:)=ZERO
     Iq=0
     DO
      Iq=Iq+1
      Call increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)
      IF (Endloop_q) exit
      W=1
      Do inb=1,size(Op%Basis%tab_basis)
       W =  W * Op%Basis%tab_basis(inb)%W(tab_iq(inb))
      END DO
      W =  W*Psi_g(Iq)**2
      DO inb=1,size(Op%Basis%tab_basis)
       Qmoy(inb) = Qmoy(inb) + Op%Basis%tab_basis(inb)%x(tab_iq(inb),1)*W
       Rho_r(iq,inb)=  W /Op%Basis%tab_basis(inb)%W(tab_iq(inb))
      END DO
     END DO
     DO inb = 1,size(Op%Basis%tab_basis)
      Write(out_unitp,*)'Qmoy(inb)',Qmoy(inb)
     END DO
     Open(1,file='Rho_r.dat',status='replace')
     DO iq = 1,Op%Basis%nq,20
     Write(1,*)iq, (Rho_r(iq,jb),jb=1,size(Op%Basis%tab_basis))
     END DO
     Close(1)
     Deallocate(Rho_r)
     Deallocate(Qmoy)
     Deallocate(tab_iq)
    END  IF
    Deallocate(EigenVal)
    Deallocate(EigenVec)

    IF (debug) THEN
     Write(out_unitp,*) 'END Diago_Op'
     flush(out_unitp)
    END IF
  END SUBROUTINE Diago_Op

  SUBROUTINE Potential(Op)
   USE Basis_m
   USE Molec_m
    Logical,          parameter         :: debug = .true.
    !Logical,         parameter         :: debug = .false.
    Logical                             :: Endloop_q
    Integer ,allocatable                :: tab_iq(:)
    Integer                             :: iq,inb
    Real (kind=Rk), allocatable         :: Q(:)
    TYPE(Op_t),  intent(inout)          :: Op
    Real (kind=Rk)                      :: V

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING Potential '
     Call Write_basis(Op%Basis)
     flush(out_unitp)
    END IF

    IF (.NOT. Basis_IS_Allocated(Op%Basis)) THEN
     Write(out_unitp,*) ' ERROR in Potential'
     Write(out_unitp,*) " the basis is not Allocated."
     STOP "ERROR  Potential: the basis is not Allocated1."
    END IF

    IF(Allocated(Op%Basis%tab_basis)) THEN
     Allocate(Tab_iq(size(Op%Basis%tab_basis)))
     Allocate(Q(size(Op%Basis%tab_basis)))
     Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
     Iq = 0
     DO
      Iq = Iq+1
      Call increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)
      IF (Endloop_q) exit
      DO inb = 1, size(Op%Basis%tab_basis)
        Q(inb) = Op%Basis%tab_basis(inb)%x(tab_iq(inb),1)
      END DO
      Call Calc_potsub(V,Q,Op%Molec)
      Op%Scalar_g(iq)= V
      Write(out_unitp,*)iq,'v=', Op%Scalar_g, 'Q',Q(:)
     END DO
     Deallocate(Tab_iq)
     Deallocate(Q)
    ELSE
     Allocate( Q(1))
     DO iq=1,Op%Basis%nq
       Q = Op%Basis%x(iq,1)
       Call Calc_potsub(V,Q,Op%Molec)
       Op%Scalar_g(iq)= V
     END DO
    END IF
    IF (debug) THEN
     Write(out_unitp,*) 'END Potential'
     flush(out_unitp)
    END IF
  END SUBROUTINE Potential

  SUBROUTINE KEO00Psi_grid(OpPsi_g,Psi_g,Op,iterm)
   USE Basis_m
   USE Molec_m
   USE UtilLib_m
    TYPE(Op_t),     intent(in),target   :: Op
    Real (kind=Rk), intent(in)          :: Psi_g(:)
    Real (kind=Rk), intent(inout),target:: OpPsi_g(:)
    Integer,        intent(in)          :: iterm
    !Logical,          parameter         :: debug = .true.
    Logical,         parameter          :: debug = .false.

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING  KEO00Psi_grid'
     Call Write_op(op)
     Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
     flush(out_unitp)
    END IF
    OpPsi_g(:) = OpPsi_g(:) +  Op%Grid(iterm)%Vec(:)* Psi_g(:)
    IF (debug) THEN
     Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     Write(out_unitp,*) 'END  KEO00Psi_grid'
     flush(out_unitp)
    END IF
  END SUBROUTINE KEO00Psi_grid


  SUBROUTINE KEOiPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq)
   USE Basis_m
   USE Molec_m
   USE UtilLib_m
    TYPE(Op_t),     intent(in),target   :: Op
    Real (kind=Rk), intent(in) ,target  :: Psi_g(:)
    Real (kind=Rk), intent(inout),target:: OpPsi_g(:)
    Integer,        intent(in)          :: inq,iterm
    Real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
    Real (kind=Rk), pointer             :: d2gg(:,:)
    Real (kind=Rk), pointer             :: d1gg(:,:)
    Real (kind=Rk), pointer             :: KEO1psi_g(:)
    Real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
    !Logical,          parameter         :: debug = .true.
    Logical,         parameter          :: debug = .false.
    Integer                             :: iq,i1,i3
    Integer                             :: Iq1,Iq2,Iq3

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING KEOiPsi_grid'
     Call Write_op(op)
     Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
     flush(out_unitp)
    END IF
    Allocate(KEO1psi_g(Op%Basis%nq))
 !!! !1) Action de d/dQi : k1_psi = d/dQi psi
    Iq1 = Product(Op%Basis%tab_basis(1:inq-1)%nq)
    Iq2 = Op%Basis%tab_basis(inq)%nq
    Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inq+1:-1)%nq)
    OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO1Psi_g
    Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g
    d1gg(1:Op%Basis%tab_basis(inq)%nq,1:Op%Basis%tab_basis(inq)%nq)=>Op%Basis%tab_basis(inq)%d1gg
    DO i3 = 1,ubound(Psi_ggg,dim=3)
    DO i1 = 1,ubound(Psi_ggg,dim=1)
      OpPsi_ggg(i1,:,i3) =  matmul(d1gg,Psi_ggg(i1,:,i3))
    END DO
    END DO
 !3) action de f1^j(Q) :  kpsi = f1^j(Q) * k1_psi
    OpPsi_g(:) = OpPsi_g(:) +  Op%Grid(iterm)%Vec(:)* KEO1Psi_g(:)
    Deallocate(KEO1psi_g)
    IF (debug) THEN
     Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     Write(out_unitp,*) 'END KEOiPsi_grid'
     flush(out_unitp)
    END IF
  END SUBROUTINE KEOiPsi_grid

 SUBROUTINE KEOijPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq1,inq2)
  USE Basis_m
  USE Molec_m
  USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   Real (kind=Rk), intent(in) ,target  :: Psi_g(:)
   Real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   Integer,        intent(in)          :: inq1,inq2,iterm
   Real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
   Real (kind=Rk), pointer             :: d2gg(:,:)
   Real (kind=Rk), pointer             :: d1gg(:,:)
   Real (kind=Rk), pointer             :: KEO1psi_g(:)
   Real (kind=Rk), pointer             :: KEO2psi_g(:)
   Real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
   !Logical,          parameter         :: debug = .true.
   Logical,         parameter          :: debug = .false.
   Integer                             :: iq,i1,i3
   Integer                             :: Iq1,Iq2,Iq3

   IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
     Call Write_op(op)
     Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
     flush(out_unitp)
   END IF
   Allocate(KEO1psi_g(Op%Basis%nq))
   Allocate(KEO2psi_g(Op%Basis%nq))

 !!!!1) Action de d/dQi : k1_psi = d/dQi psi
   Iq1 = Product(Op%Basis%tab_basis(1:inq1-1)%nq)
   Iq2 = Op%Basis%tab_basis(inq1)%nq
   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inq1+1:-1)%nq)
   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO1Psi_g
   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g
   d1gg(1:Op%Basis%tab_basis(inq1)%nq,1:Op%Basis%tab_basis(inq1)%nq) => Op%Basis%tab_basis(inq1)%d1gg
   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)
     OpPsi_ggg(i1,:,i3) = matmul(d1gg,Psi_ggg(i1,:,i3))
   END DO
   END DO
!!!!2) Action de d/dQj :   k2_psi = d/dQj k1_psi
   Iq1 = Product(Op%Basis%tab_basis(1:inq2-1)%nq)
   Iq2 = Op%Basis%tab_basis(inq2)%nq
   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inq2+1:-1)%nq)
   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO2Psi_g
   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => KEO1Psi_g
   d1gg(1:Op%Basis%tab_basis(inq2)%nq,1:Op%Basis%tab_basis(inq2)%nq)=>Op%Basis%tab_basis(inq2)%d1gg
   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)
     OpPsi_ggg(i1,:,i3) =  matmul(d1gg,Psi_ggg(i1,:,i3))
   END DO
   END DO
!3) action de f2^ij(Q) :  kpsi = f2^ij(Q) * k2_psi
   OpPsi_g(:) = OpPsi_g(:) +   Op%Grid(iterm)%Vec(:)*KEO2Psi_g(:)
   Deallocate(KEO1psi_g)
   Deallocate(KEO2psi_g)
   IF (debug) THEN
     Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     Write(out_unitp,*) 'END KEOijPsi_grid'
     flush(out_unitp)
   END IF
 END SUBROUTINE KEOijPsi_grid


 SUBROUTINE KEOiiPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq)
  USE Basis_m
  USE Molec_m
  USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   Real (kind=Rk), intent(in) ,target  :: Psi_g(:)
   Real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   Integer,        intent(in)          :: inq,iterm
   Real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
   Real (kind=Rk), pointer             :: Psi_ggg1(:,:,:)
   Real (kind=Rk), pointer             :: d2gg(:,:)
   Real (kind=Rk), pointer             :: d1gg(:,:)
   Real (kind=Rk), pointer             :: KEO1psi_g(:)
   Real (kind=Rk), pointer             :: KEO2psi_g(:),KEO21Psi_g(:)
   Real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
   Real (kind=Rk), pointer             :: OpPsi_ggg1(:,:,:)
   !Logical,          parameter         :: debug = .true.
   Logical,         parameter          :: debug = .false.
   Integer                             :: iq,i1,i3
   Integer                             :: Iq1,Iq2,Iq3

   IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
     Call Write_op(op)
     Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
     flush(out_unitp)
   END IF
   Allocate(KEO1psi_g(Op%Basis%nq))
   Allocate(KEO2psi_g(Op%Basis%nq))
   Allocate(KEO21psi_g(Op%Basis%nq))
   !!!!1) Action de d/dQi : k1_psi = d/dQi psi
   Iq1 = Product(Op%Basis%tab_basis(1:inq-1)%nq)
   Iq2 = Op%Basis%tab_basis(inq)%nq
   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inq+1:-1)%nq)
   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO2Psi_g
   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g
   d1gg(1:Op%Basis%tab_basis(inq)%nq,1:Op%Basis%tab_basis(inq)%nq)=>Op%Basis%tab_basis(inq)%d1gg

   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)
     OpPsi_ggg(i1,:,i3) = matmul(d1gg,  matmul(d1gg,Psi_ggg(i1,:,i3)))
   END DO
   END DO
	 !3) action de f2^ij(Q) :  kpsi = f2^ij(Q) * k2_psi
   OpPsi_g(:) = OpPsi_g(:) +  Op%Grid(iterm)%Vec(:)*KEO2Psi_g(:)
   Deallocate(KEO1psi_g)
   Deallocate(KEO2psi_g)
   Deallocate(KEO21psi_g)
   IF (debug) THEN
	   Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
	   Write(out_unitp,*) 'END KEOiiPsi_grid'
	   flush(out_unitp)
   END IF
 END SUBROUTINE KEOiiPsi_grid


 SUBROUTINE KEOii1Psi_grid(OpPsi_g,Psi_g,Op,iterm,inq)
  USE Basis_m
  USE Molec_m
  USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   Real (kind=Rk), intent(in) ,target  :: Psi_g(:)
   Real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   Integer,        intent(in)          :: inq,iterm
   Real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
   Real (kind=Rk), pointer             :: d2gg(:,:)
   Real (kind=Rk), pointer             :: d1gg(:,:)
   Real (kind=Rk), pointer             :: KEO2psi_g(:)
   Real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
   !Logical,          parameter         :: debug = .true.
   Logical,         parameter          :: debug = .false.
   Integer                             :: iq,i1,i3
   Integer                             :: Iq1,Iq2,Iq3

   IF (debug) THEN
    Write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
    Call Write_op(op)
    Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
   END IF
   Allocate(KEO2psi_g(Op%Basis%nq))
   Iq1 = Product(Op%Basis%tab_basis(1:inq-1)%nq)
   Iq2 = Op%Basis%tab_basis(inq)%nq
   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inq+1:-1)%nq)
   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO2Psi_g
   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g
   d2gg(1:Op%Basis%tab_basis(inq)%nq,1:Op%Basis%tab_basis(inq)%nq)=>Op%Basis%tab_basis(inq)%d2gg

   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)
     OpPsi_ggg(i1,:,i3) =   matmul(d2gg,Psi_ggg(i1,:,i3))
   END DO
   END DO
   OpPsi_g(:) = OpPsi_g(:) +  Op%Grid(iterm)%Vec(:)*KEO2Psi_g(:)
   Deallocate(KEO2psi_g)
   IF (debug) THEN
     Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     Write(out_unitp,*) 'END KEOiiPsi_grid'
     flush(out_unitp)
   END IF
 END SUBROUTINE KEOii1Psi_grid

 SUBROUTINE OpPsi_gridcnD(OpPsi_g,Psi_g,Op)
  USE Basis_m
  USE Molec_m
  USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   Real (kind=Rk), intent(in) ,target  :: Psi_g(:)
   Real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   Real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
   Real (kind=Rk), pointer             :: d2gg(:,:)
   Real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
   !Logical,          parameter         :: debug = .true.
   Logical,         parameter          :: debug = .false.
   Integer                             :: iq,i1,i3,inq1,inq2
   Integer                             :: Iq1,Iq2,Iq3,iterm

   IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING OpPsi_gridcnD'
     Call Write_op(op)
     Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
     flush(out_unitp)
   END IF

   OpPsi_g  = ZERO
   DO iterm = 1,size(Op%Grid)
     inq1  = Op%Grid(iterm)%DerivIndex(1)
     inq2  = Op%Grid(iterm)%DerivIndex(2)
     IF (inq1 > 0 .and. inq2 > 0 .and. inq1 /= inq2) THEN ! terme croisé
       Call KEOijPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq1,inq2)
     ELSE IF (inq1 > 0 .and. inq2 > 0 .and. inq1 == inq2) THEN ! terme non croisé
       Call KEOii1Psi_grid(OpPsi_g,Psi_g,Op,iterm,inq1)
     ELSE IF (inq1 > 0 .and. inq2 == 0) THEN
       Call KEOiPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq1) ! dérivé premiere
     ELSE IF (inq1 == 0 .and. inq2 == 0) THEN
       Call KEO00Psi_grid(OpPsi_g,Psi_g,Op,iterm)     ! les scalaires
     END IF
   END DO

   IF (debug) THEN
     Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     Write(out_unitp,*) 'END OpPsi_gridcnD'
     flush(out_unitp)
   END IF
 END SUBROUTINE OpPsi_gridcnD


 SUBROUTINE OpPsi_gridnD(OpPsi_g,Psi_g,Op)
  USE Basis_m
  USE Molec_m
  USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   Real (kind=Rk), intent(in), target  :: Psi_g(:)
   Real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   Real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
   Real (kind=Rk), pointer             :: d2gg(:,:)
   Real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
   !Logical,          parameter         :: debug = .true.
   Logical,         parameter          :: debug = .false.
   Integer                             :: iq,i1,i3,inb
   Integer                             :: Iq1,Iq2,Iq3

   IF (debug) THEN
    Write(out_unitp,*) 'BEGINNING OpPsi_gridnD'
    Call Write_op(op)
    Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
   END IF

   OpPsi_g(:) = Op%Scalar_g(:)*Psi_g(:)

   DO inb = 1,size(Op%Basis%tab_basis)
    Iq1 = Product(Op%Basis%tab_basis(1:inb-1)%nq)
    Iq2 = Op%Basis%tab_basis(inb)%nq
    Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inb+1:-1)%nq)
    OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> OpPsi_g
    Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g
    d2gg(1:Op%Basis%tab_basis(inb)%nq,1:Op%Basis%tab_basis(inb)%nq)=>Op%Basis%tab_basis(inb)%d2gg
    DO i3 = 1,ubound(Psi_ggg,dim=3)
    DO i1 = 1,ubound(Psi_ggg,dim=1)
      OpPsi_ggg(i1,:,i3) = OpPsi_ggg(i1,:,i3) -HALF/mass* matmul(d2gg,Psi_ggg(i1,:,i3))
    END DO
    END DO
   END DO

   IF (debug) THEN
    Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
    Write(out_unitp,*) 'END OpPsi_gridnD'
    flush(out_unitp)
   END IF

 END SUBROUTINE OpPsi_gridnD

 SUBROUTINE OpPsi_grid(OpPsi_g,Psi_g,Op)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
!  TYPE (Molec_t),  intent(in)         :: Molec
  TYPE(Op_t) , intent(inout)          :: Op
  Real (kind=Rk), intent(in) ,target  :: Psi_g(:)
  Real (kind=Rk), intent(inout)       :: OpPsi_g(:)
  !Logical,          parameter         :: debug = .true.
  Logical,         parameter          :: debug  = .false.
  Logical,         parameter          :: curv   = .true.
  Logical,         parameter          :: simple = .false.

  IF (debug) THEN
    Write(out_unitp,*) 'BEGINNING OpPsi_grid'
    Call Write_op(op)
    Call Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
  END IF

  IF (Allocated(Op%Basis%tab_basis)) THEN
    IF (curv) THEN
      Call OpPsi_gridcnD(OpPsi_g,Psi_g,Op)
    ELSE IF (Simple) THEN
      IF (.NOT.Allocated(Op%Scalar_g)) THEN
        Allocate(Op%Scalar_g(Op%Basis%nq))
        Call Potential(Op)
      END IF
      Call OpPsi_gridnD(OpPsi_g,Psi_g,Op)
    END IF
  ELSE IF (Basis_IS_Allocated(Op%Basis)) THEN
   OpPsi_g(:) = Op%Scalar_g(:) * Psi_g(:)
   OpPsi_g(:) = OpPsi_g(:) -HALF/mass * matmul(Op%Basis%d2gg(:,:,1,1),Psi_g(:))
  END IF
  IF (debug) THEN
    Call Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
    Write(out_unitp,*) 'END OpPsi_grid'
    flush(out_unitp)
  END IF
 END SUBROUTINE OpPsi_grid

 SUBROUTINE calc_OpPsi(Op,Psi,OpPsi)
 USE psi_m, ONLY : psi_t
  TYPE(Op_t),  intent(in)     :: Op
  TYPE(psi_t), intent(in)     :: Psi
  TYPE(psi_t), intent(inout)  :: OpPsi

  IF (.NOT. Allocated(Op%RMat)) THEN
    STOP 'ERROR in calc_OpPsi: Op is not initialized!'
  END IF

  IF (Allocated(Psi%RVec)) THEN
    OpPsi%RVec = matmul(Op%RMat,Psi%RVec)
  ELSE IF (Allocated(Psi%CVec)) THEN
    OpPsi%CVec = matmul(Op%RMat,Psi%CVec)
  ELSE
    STOP 'ERROR in calc_OpPsi: Psi is not initialized!'
  END IF

 END SUBROUTINE calc_OpPsi
END MODULE Op_m
