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
module Op_m
  USE NumParameters_m
  USE diago_m
  USE NDindex_m
  USE Molec_m
  USE Basis_m, only : Basis_t

  implicit none
  private

  TYPE :: Grid_t
    REAL(kind=Rk),    allocatable :: Vec(:)
    integer                       :: DerivIndex(2)
  END TYPE Grid_t

  TYPE :: Op_t
    TYPE (Grid_t) ,    allocatable :: Grid(:)
    TYPE (Basis_t),    pointer     :: Basis => Null()
    TYPE (Molec_t),    pointer     :: Molec => Null()
    real (kind=Rk),    allocatable :: Scalar_g(:) ! The scalar part of the operator (The  for Hamiltonian)
    real (kind=Rk),    allocatable :: RMat(:,:) ! The matrix of the Hamiltonian operator

  END TYPE Op_t



  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi,Diago_Op,Make_Mat_Op
  public :: OpPsi_gridnD,Grid_t

contains
  SUBROUTINE alloc_Op(Op,nb)
    TYPE(Op_t),  intent(inout) :: Op
    integer,     intent(in)    :: nb

    CALL dealloc_Op(Op)

    IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

    allocate(Op%RMat(nb,nb))

  END SUBROUTINE alloc_Op

  SUBROUTINE dealloc_Op(Op)
    TYPE(Op_t), intent(inout) :: Op

    nullify(Op%Basis)

    IF (allocated(Op%RMat)) THEN
      deallocate(Op%RMat)
    END IF

    IF (allocated(Op%Scalar_g)) THEN
      deallocate(Op%Scalar_g)
    END IF

  END SUBROUTINE dealloc_Op

  SUBROUTINE write_Op(Op)
  USE UtilLib_m, only : Write_RMat

   TYPE(Op_t), intent(in) :: Op

   IF (associated(Op%Basis)) THEN
    write(out_unitp,*) ' The basis is linked to Op.'
   ELSE
    write(out_unitp,*) ' The basis is NOT linked to Op.'
   END IF

   IF (allocated(Op%RMat)) THEN
     write(out_unitp,*) 'Writing Op (real):'
     write(out_unitp,*)
     CALL Write_RMat(Op%RMat,out_unitp,5,name_info='Op%Rmat')
     write(out_unitp,*) 'END Writing Op'
   END IF

 END SUBROUTINE write_Op


 SUBROUTINE Set_Op(Op,Basis,Molec)
 USE Basis_m
 USE Molec_m
 USE NDindex_m

  TYPE(Op_t),     intent(inout)       :: Op
  TYPE (Basis_t), intent(in),  target :: Basis
  TYPE (Molec_t), intent(in),  target :: Molec
  !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: ib,iq
  logical                             :: Simple = .false.
!----- for debuging --------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING Set_Op'
    call write_basis(Basis)
    flush(out_unitp)
  END IF
!----------------------end --------------------------------------------
  Op%Basis => Basis  !Initialization of Op%Basis by the values ​​of Basis
  Op%Molec => Molec

  IF (Simple) THEN
   allocate(Op%Grid(3))

   allocate(Op%Grid(1)%Vec(Op%Basis%nq))
   allocate(Op%Grid(2)%Vec(Op%Basis%nq))
   allocate(Op%Grid(3)%Vec(Op%Basis%nq))

   Op%Grid(1)%DerivIndex(1)=1
   Op%Grid(1)%DerivIndex(2)=1
   Op%Grid(2)%DerivIndex(1)=2
   Op%Grid(2)%DerivIndex(2)=2
   Op%Grid(3)%DerivIndex(1)=3
   Op%Grid(3)%DerivIndex(2)=3


   Op%Grid(1)%Vec(:)=-HALF/mass
   Op%Grid(2)%Vec(:)=-HALF/mass
   Op%Grid(3)%Vec(:)=-HALF/mass
  ELSE
   CALL Set_grid_Op(Op)
  END IF
  IF (debug) THEN
    write(out_unitp,*) 'END Set_Op'
    flush(out_unitp)
  END IF
 END SUBROUTINE Set_Op

 SUBROUTINE Set_grid_Op(Op)
 USE Basis_m
 USE NDindex_m

  TYPE(Op_t),     intent(inout)       :: Op
  real(kind=Rk)                       :: F,F2(3,3),F1(3),Vep,V
  real(kind=Rk), allocatable          :: Q(:)
  !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: ib,iq,inb
  logical                             :: Endloop_q
  logical                             :: Endloop_b
  integer,         allocatable        :: tab_iq(:)


  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING Set_grid_Op'
    call write_basis(Op%Basis)
    flush(out_unitp)
  END IF
  allocate(Q(size(Op%Basis%tab_basis)))
  allocate(Op%Grid(14))
  allocate(Op%Grid(1)%Vec(Op%Basis%nq))
  allocate(Op%Grid(2)%Vec(Op%Basis%nq))
  allocate(Op%Grid(3)%Vec(Op%Basis%nq))
  allocate(Op%Grid(4)%Vec(Op%Basis%nq))
  allocate(Op%Grid(5)%Vec(Op%Basis%nq))
  allocate(Op%Grid(6)%Vec(Op%Basis%nq))
  allocate(Op%Grid(7)%Vec(Op%Basis%nq))
  allocate(Op%Grid(8)%Vec(Op%Basis%nq))
  allocate(Op%Grid(9)%Vec(Op%Basis%nq))
  allocate(Op%Grid(10)%Vec(Op%Basis%nq))
  allocate(Op%Grid(11)%Vec(Op%Basis%nq))
  allocate(Op%Grid(12)%Vec(Op%Basis%nq))
  allocate(Op%Grid(13)%Vec(Op%Basis%nq))
  allocate(Op%Grid(14)%Vec(Op%Basis%nq))

  Op%Grid(1)%DerivIndex(1)  = 0
  Op%Grid(1)%DerivIndex(2)  = 0
  Op%Grid(2)%DerivIndex(1)  = 1
  Op%Grid(2)%DerivIndex(2)  = 1
  Op%Grid(3)%DerivIndex(1)  = 1
  Op%Grid(3)%DerivIndex(2)  = 2
  Op%Grid(4)%DerivIndex(1)  = 1
  Op%Grid(4)%DerivIndex(2)  = 3
  Op%Grid(5)%DerivIndex(1)  = 2
  Op%Grid(5)%DerivIndex(2)  = 1
  Op%Grid(6)%DerivIndex(1)  = 2
  Op%Grid(6)%DerivIndex(2)  = 2
  Op%Grid(7)%DerivIndex(1)  = 2
  Op%Grid(7)%DerivIndex(2)  = 3
  Op%Grid(8)%DerivIndex(1)  = 3
  Op%Grid(8)%DerivIndex(2)  = 1
  Op%Grid(9)%DerivIndex(1)  = 3
  Op%Grid(9)%DerivIndex(2)  = 2
  Op%Grid(10)%DerivIndex(1) = 3
  Op%Grid(10)%DerivIndex(2) = 3
  Op%Grid(11)%DerivIndex(1) = 1
  Op%Grid(11)%DerivIndex(2) = 0
  Op%Grid(12)%DerivIndex(1) = 2
  Op%Grid(12)%DerivIndex(2) = 0
  Op%Grid(13)%DerivIndex(1) = 3
  Op%Grid(13)%DerivIndex(2) = 0
  Op%Grid(14)%DerivIndex(1) = 0
  Op%Grid(14)%DerivIndex(2) = 0


  Allocate(Tab_iq(size(Op%Basis%tab_basis)))
  Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
  Iq=0
  DO
   Iq=Iq+1
   CALL increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)

   IF (Endloop_q) exit
   Do inb = 1,size(Op%Basis%tab_basis)
     Q(inb) =  Op%Basis%tab_basis(inb)%x(tab_iq(inb))
   END DO
   Call Tana_F2_F1_Vep(F2,F1,Vep,Q)
   CALL Calc_potsub(V,Q,Op%Molec)
   Op%Grid(1)%Vec(Iq)  = V
   Op%Grid(2)%Vec(Iq)  = F2(1,1)
   Op%Grid(3)%Vec(Iq)  = F2(1,2)
   Op%Grid(4)%Vec(Iq)  = F2(1,3)
   Op%Grid(5)%Vec(Iq)  = F2(2,1)
   Op%Grid(6)%Vec(Iq)  = F2(2,2)
   Op%Grid(7)%Vec(Iq)  = F2(2,3)
   Op%Grid(8)%Vec(Iq)  = F2(3,1)
   Op%Grid(9)%Vec(Iq)  = F2(3,2)
   Op%Grid(10)%Vec(Iq) = F2(3,3)
   Op%Grid(11)%Vec(Iq) = F1(1)
   Op%Grid(12)%Vec(Iq) = F1(2)
   Op%Grid(13)%Vec(Iq) = F1(3)
   Op%Grid(14)%Vec(Iq) = Vep
  END DO

  Deallocate(Tab_iq)

  IF (debug) THEN
    write(out_unitp,*) 'END Set_drid_Op'
    flush(out_unitp)
  END IF

 END SUBROUTINE Set_grid_op


 SUBROUTINE Make_Mat_OP(Op)
 USE Molec_m
 USE Basis_m

   TYPE (Op_t),     intent(inout)      :: Op
   !logical,          parameter         :: debug = .true.
   logical,         parameter          :: debug = .false.
   integer                             :: ib,iq,jb
   real (kind=Rk), allocatable         :: Psi_g(:)
   real (kind=Rk), allocatable         :: Psi_b(:)!,Psi_g(:)
   real (kind=Rk), allocatable         :: OpPsi_g(:)

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING Make_Mat_OP'
     call write_Op(Op)
     call write_basis(Op%Basis)
     flush(out_unitp)
   END IF

   allocate(Op%RMat(Op%Basis%nb,Op%Basis%nb))
   allocate(Psi_b(Op%Basis%nb))
   allocate(Psi_g(Op%Basis%nq))
   allocate(OpPsi_g(Op%Basis%nq))

   DO ib=1,Op%Basis%nb
     Psi_b(:)  = ZERO
     Psi_b(ib) = ONE

     CALL BasisTOGrid_Basis(Psi_g, Psi_b,Op%Basis)
     CALL OpPsi_grid(OpPsi_g,Psi_g,Op)
     CALL GridTOBasis_Basis(Op%RMat(:,ib), OpPsi_g,Op%Basis)
   END DO

   Deallocate(Psi_b)
   Deallocate(Psi_g)
   Deallocate(OpPsi_g)

   IF (debug) THEN
     write(out_unitp,*) 'END Make_Mat_OP'
     flush(out_unitp)
   END IF

  END SUBROUTINE Make_Mat_OP

  SUBROUTINE Diago_Op(Op)
  USE Molec_m
  USE Basis_m
   TYPE(Op_t),     intent(inout)       :: Op
   !logical,          parameter         :: debug = .true.
   logical,         parameter          ::debug = .false.
   integer                             :: jb,inb
   integer                             :: ib,iq
   logical                             :: Endloop_q
   integer,         allocatable        :: tab_iq(:),QW(:,:),Qmoy(:),Rho_r(:,:)
   real (kind=Rk), allocatable         :: EigenVal(:),EigenVec(:,:),Psi_g(:),DifEigen(:)

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING Diago_Op'
     flush(out_unitp)
   END IF
   allocate(Psi_g(Op%Basis%nq))
   allocate(DifEigen(Op%Basis%nb))
   allocate(EigenVal(Op%Basis%nb))
   allocate(EigenVec(Op%Basis%nb,Op%Basis%nb))

   CALL  diagonalization(Op%RMat,EigenVal,EigenVec,Op%Basis%nb)

   Write(out_unitp,*)
   Write(out_unitp,*)
   Write(out_unitp,*) 'eigenvalues = '
   write(out_unitp,*) 'n','EigenVal','E_n-E_1'
   DO ib=1,5!Op%Basis%nb
     DifEigen(ib)=(EigenVal(ib)-Op%Molec%V0)*219474.631443_RK-(EigenVal(1)-Op%Molec%V0)*219474.631443_RK
     write(out_unitp,*) ib,(EigenVal(ib)-Op%Molec%V0)*219474.631443_RK,'DifEigen(ib)',DifEigen(ib)
   END DO

   Write(out_unitp,*)
   Write(out_unitp,*)

   CALL BasisTOGrid_Basis(Psi_g, EigenVec(:,1),Op%Basis)

   Allocate(Qmoy(size(Op%Basis%tab_basis)))
   Allocate(Rho_r(Op%Basis%nq,size(Op%Basis%tab_basis)))
   Allocate(Qw(Op%Basis%nq,size(Op%Basis%tab_basis)))
   Allocate(Tab_iq(size(Op%Basis%tab_basis)))
   Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
    Iq=0
    DO
    Iq=Iq+1
    CALL increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)
    IF (Endloop_q) exit
      Do inb=1,size(Op%Basis%tab_basis)
        Qw(Iq,inb) =  Op%Basis%tab_basis(inb)%x(tab_iq(inb))*&
                              Op%Basis%tab_basis(inb)%W(tab_iq(inb))
      END DO
      Rho_r(Iq,1) = Psi_g(Iq)*Psi_g(Iq) *Op%Basis%tab_basis(2)%W(tab_iq(2))&
                                *Op%Basis%tab_basis(3)%W(tab_iq(3))
      Rho_r(Iq,2) = Psi_g(Iq)*Psi_g(Iq) *Op%Basis%tab_basis(1)%W(tab_iq(1))&
                                *Op%Basis%tab_basis(3)%W(tab_iq(3))
      Rho_r(Iq,3) = Psi_g(Iq)*Psi_g(Iq) *Op%Basis%tab_basis(1)%W(tab_iq(1))&
                                *Op%Basis%tab_basis(2)%W(tab_iq(2))
    END DO
    DO inb = 1,size(Op%Basis%tab_basis)
    Qmoy(inb) = dot_product(Psi_g(:)*Psi_g(:),Qw(:,inb))
    Write(out_unitp,*)'Qmoy(inb)',Qmoy(inb)
   END DO

   Open(1,file='Rho_r.dat',status='replace')
   Do inb = 1,Op%Basis%nq  
      write(1,*)(Rho_r(inb,jb),jb=1,size(Op%Basis%tab_basis))
   END DO
   Close(1)


   Deallocate(Rho_r)
   Deallocate(QW)
   Deallocate(Qmoy)
   Deallocate(tab_iq)
   Deallocate(EigenVal)
   Deallocate(EigenVec)

   IF (debug) THEN
     write(out_unitp,*) 'END Diago_Op'
     flush(out_unitp)
   END IF

  END SUBROUTINE Diago_Op

  SUBROUTINE Potential(Op)
  USE Basis_m
  USE Molec_m
   logical,          parameter         :: debug = .true.
   !logical,         parameter          :: debug = .false.
   logical                             :: Endloop_q
   integer ,allocatable                :: tab_iq(:)
   integer                             :: iq,inb
   real (kind=Rk), allocatable         :: Q(:)
   TYPE(Op_t),  intent(inout)          :: Op
   real (kind=Rk)                      :: V


   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING Potential '
     CALL Write_basis(Op%Basis)
     flush(out_unitp)
   END IF

   IF (.NOT. Basis_IS_allocated(Op%Basis)) THEN
     write(out_unitp,*) ' ERROR in Potential'
     write(out_unitp,*) " the basis is not allocated."
     STOP "ERROR  Potential: the basis is not allocated1."
   END IF

   IF(allocated(Op%Basis%tab_basis)) THEN
     allocate(Tab_iq(size(Op%Basis%tab_basis)))
     allocate(Q(size(Op%Basis%tab_basis)))
     Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
     Iq = 0
     DO
      Iq = Iq+1
      CALL increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)
      IF (Endloop_q) exit
      DO inb = 1, size(Op%Basis%tab_basis)
        Q(inb) = Op%Basis%tab_basis(inb)%x(tab_iq(inb))
      END DO
      CALL Calc_potsub(V,Q,Op%Molec)
      Op%Scalar_g(iq)= V
      write(out_unitp,*)iq,'v=', Op%Scalar_g, 'Q',Q(:)
      !Op%Scalar_g(iq)=Calc_pot(Q)
     END DO
     Deallocate(Tab_iq)
     Deallocate(Q)
   ELSE
     allocate( Q(1))
     DO iq=1,Op%Basis%nq
       Q = Op%Basis%x(iq)
       CALL Calc_potsub(V,Q,Op%Molec)
       Op%Scalar_g(iq)= V
       !Op%Scalar_g(iq) = Calc_pot(Q)
     END DO
   END IF

   IF (debug) THEN

     write(out_unitp,*) 'END Potential'
     flush(out_unitp)
   END IF
 END SUBROUTINE Potential

 SUBROUTINE KEO00Psi_grid(OpPsi_g,Psi_g,Op,iterm)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   real (kind=Rk), intent(in)          :: Psi_g(:)
   real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   integer,        intent(in)          :: iterm
  !logical,          parameter         :: debug = .true.
   logical,         parameter          :: debug = .false.


   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING  KEO00Psi_grid'
     CALL Write_op(op)
     CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
     flush(out_unitp)
   END IF


   OpPsi_g(:) = OpPsi_g(:) +  Op%Grid(iterm)%Vec(:)* Psi_g(:)


   IF (debug) THEN
     CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     write(out_unitp,*) 'END  KEO00Psi_grid'
     flush(out_unitp)
   END IF

 END SUBROUTINE KEO00Psi_grid


 SUBROUTINE KEOiPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
  TYPE(Op_t),     intent(in),target   :: Op
  real (kind=Rk), intent(in) ,target  :: Psi_g(:)
  real (kind=Rk), intent(inout),target:: OpPsi_g(:)
  integer,        intent(in)          :: inq,iterm
  real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
  real (kind=Rk), pointer             :: d2gg(:,:)
  real (kind=Rk), pointer             :: d1gg(:,:)
  real (kind=Rk), pointer             :: KEO1psi_g(:)
  real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
 !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: iq,i1,i3
  integer                             :: Iq1,Iq2,Iq3

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING KEOiPsi_grid'
    CALL Write_op(op)
    CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
  END IF

 !!!!1) Action de d/dQi : k1_psi = d/dQi psi

  Iq1 = Product(Op%Basis%tab_basis(1:inq-1)%nq)

  Iq2 = Op%Basis%tab_basis(inq)%nq

  Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inq+1:-1)%nq)

  OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO1Psi_g

  Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g

  d1gg(1:Op%Basis%tab_basis(inq)%nq,1:Op%Basis%tab_basis(inq)%nq)=>Op%Basis%tab_basis(inq)%d1gg

  DO i3=1,ubound(Psi_ggg,dim=3)
  DO i1=1,ubound(Psi_ggg,dim=1)
     OpPsi_ggg(i1,:,i3) =  matmul(d1gg,Psi_ggg(i1,:,i3))
  END DO
  END DO

 !3) action de f1^j(Q) :  kpsi = f1^j(Q) * k1_psi

  OpPsi_g(:) = OpPsi_g(:) +  Op%Grid(iterm)%Vec(:)* KEO1Psi_g(:)

  IF (debug) THEN
    CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
    write(out_unitp,*) 'END KEOiPsi_grid'
    flush(out_unitp)
  END IF

 END SUBROUTINE KEOiPsi_grid

 SUBROUTINE KEOijPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq1,inq2)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
  TYPE(Op_t),     intent(in),target   :: Op
  real (kind=Rk), intent(in) ,target  :: Psi_g(:)
  real (kind=Rk), intent(inout),target:: OpPsi_g(:)
  integer,        intent(in)          :: inq1,inq2,iterm
  real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
  real (kind=Rk), pointer             :: d2gg(:,:)
  real (kind=Rk), pointer             :: d1gg(:,:)
  real (kind=Rk), pointer             :: KEO1psi_g(:)
  real (kind=Rk), pointer             :: KEO2psi_g(:)
  real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
 !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: iq,i1,i3
  integer                             :: Iq1,Iq2,Iq3

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
    CALL Write_op(op)
    CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
  END IF

 !!!!1) Action de d/dQi : k1_psi = d/dQi psi

  Iq1 = Product(Op%Basis%tab_basis(1:inq1-1)%nq)

  Iq2 = Op%Basis%tab_basis(inq1)%nq

  Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inq1+1:-1)%nq)

  OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO1Psi_g

  Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g

  d1gg(1:Op%Basis%tab_basis(inq1)%nq,1:Op%Basis%tab_basis(inq1)%nq)=>Op%Basis%tab_basis(inq1)%d1gg

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


  IF (debug) THEN
     CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     write(out_unitp,*) 'END KEOijPsi_grid'
     flush(out_unitp)
  END IF

 END SUBROUTINE KEOijPsi_grid


 SUBROUTINE KEOiiPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
  TYPE(Op_t),     intent(in),target   :: Op
  real (kind=Rk), intent(in) ,target  :: Psi_g(:)
  real (kind=Rk), intent(inout),target:: OpPsi_g(:)
  integer,        intent(in)          :: inq,iterm
  real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
  real (kind=Rk), pointer             :: Psi_ggg1(:,:,:)
  real (kind=Rk), pointer             :: d2gg(:,:)
  real (kind=Rk), pointer             :: d1gg(:,:)
  real (kind=Rk), pointer             :: KEO1psi_g(:)
  real (kind=Rk), pointer             :: KEO2psi_g(:),KEO21Psi_g(:)
  real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
  real (kind=Rk), pointer             :: OpPsi_ggg1(:,:,:)
 !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: iq,i1,i3
  integer                             :: Iq1,Iq2,Iq3

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
    CALL Write_op(op)
    CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
  END IF

  allocate(KEO1psi_g(Op%Basis%nq))
  allocate(KEO2psi_g(Op%Basis%nq))
  allocate(KEO21psi_g(Op%Basis%nq))

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


     IF (debug) THEN
	  CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
	  write(out_unitp,*) 'END KEOiiPsi_grid'
	  flush(out_unitp)
     END IF

  END SUBROUTINE KEOiiPsi_grid



 SUBROUTINE KEOii1Psi_grid(OpPsi_g,Psi_g,Op,iterm,inq)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
  TYPE(Op_t),     intent(in),target   :: Op
  real (kind=Rk), intent(in) ,target  :: Psi_g(:)
  real (kind=Rk), intent(inout),target:: OpPsi_g(:)
  integer,        intent(in)          :: inq,iterm
  real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
  real (kind=Rk), pointer             :: d2gg(:,:)
  real (kind=Rk), pointer             :: d1gg(:,:)
  real (kind=Rk), pointer             :: KEO2psi_g(:)
  real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
  !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: iq,i1,i3
  integer                             :: Iq1,Iq2,Iq3

  IF (debug) THEN
   write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
   CALL Write_op(op)
   CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
   flush(out_unitp)
  END IF

  allocate(KEO2psi_g(Op%Basis%nq))

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


  IF (debug) THEN
   CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
   write(out_unitp,*) 'END KEOiiPsi_grid'
   flush(out_unitp)
  END IF

 END SUBROUTINE KEOii1Psi_grid

 SUBROUTINE OpPsi_gridcnD(OpPsi_g,Psi_g,Op)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
  TYPE(Op_t),     intent(in),target   :: Op
  real (kind=Rk), intent(in) ,target  :: Psi_g(:)
  real (kind=Rk), intent(inout),target:: OpPsi_g(:)
  real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
  real (kind=Rk), pointer             :: d2gg(:,:)
  real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
  !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: iq,i1,i3,inq1,inq2
  integer                             :: Iq1,Iq2,Iq3,iterm


  IF (debug) THEN
   write(out_unitp,*) 'BEGINNING OpPsi_gridcnD'
   CALL Write_op(op)
   CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
   flush(out_unitp)
  END IF

  !OpPsi_g(:) = Op%Scalar_g(:)*Psi_g(:)
  OpPsi_g = 0._Rk
  DO iterm=1,size(Op%Grid)
   inq1 = Op%Grid(iterm)%DerivIndex(1)
   inq2 = Op%Grid(iterm)%DerivIndex(2)

   IF (inq1 > 0 .and. inq2 > 0 .and. inq1 /= inq2) THEN ! terme croisé
    CALL KEOijPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq1,inq2)
   ELSE IF (inq1 > 0 .and. inq2 > 0 .and. inq1 == inq2) THEN ! terme non croisé
    CALL KEOiiPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq1)
   ELSE IF (inq1 > 0 .and. inq2 == 0) THEN
    CALL KEOiPsi_grid(OpPsi_g,Psi_g,Op,iterm,inq1) ! dérivé premiere
   ELSE IF (inq1 == 0 .and. inq2 == 0) THEN
    CALL KEO00Psi_grid(OpPsi_g,Psi_g,Op,iterm)
   END IF
  END DO

  IF (debug) THEN
   CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
   write(out_unitp,*) 'END OpPsi_gridcnD'
   flush(out_unitp)
  END IF

 END SUBROUTINE OpPsi_gridcnD


 SUBROUTINE OpPsi_gridnD(OpPsi_g,Psi_g,Op)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
  TYPE(Op_t),     intent(in),target   :: Op
  real (kind=Rk), intent(in), target  :: Psi_g(:)
  real (kind=Rk), intent(inout),target:: OpPsi_g(:)
  real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
  real (kind=Rk), pointer             :: d2gg(:,:)
  real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
  !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  integer                             :: iq,i1,i3,inb
  integer                             :: Iq1,Iq2,Iq3

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING OpPsi_gridnD'
    CALL Write_op(op)
    CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
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
    DO i3=1,ubound(Psi_ggg,dim=3)
    DO i1=1,ubound(Psi_ggg,dim=1)
      OpPsi_ggg(i1,:,i3) = OpPsi_ggg(i1,:,i3) -HALF/mass* matmul(d2gg,Psi_ggg(i1,:,i3))
    END DO
    END DO
  END DO

  IF (debug) THEN
    CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
    write(out_unitp,*) 'END OpPsi_gridnD'
    flush(out_unitp)
  END IF

 END SUBROUTINE OpPsi_gridnD


 SUBROUTINE OpPsi_grid(OpPsi_g,Psi_g,Op)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
!  TYPE (Molec_t),  intent(in)         :: Molec
  TYPE(Op_t) , intent(inout)          :: Op
  real (kind=Rk), intent(in) ,target  :: Psi_g(:)
  real (kind=Rk), intent(inout)       :: OpPsi_g(:)
  !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  logical,         parameter          :: curv = .true.
  logical,         parameter          :: simple = .false.

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING OpPsi_grid'
    CALL Write_op(op)
    CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
  END IF

  IF (allocated(Op%Basis%tab_basis)) THEN
    IF (curv) THEN
      CALL OpPsi_gridcnD(OpPsi_g,Psi_g,Op)
    ELSE IF (Simple) THEN
      IF (.NOT.allocated(Op%Scalar_g)) THEN
        ALLOCATE(Op%Scalar_g(Op%Basis%nq))
        CALL Potential(Op)
      END IF
      CALL OpPsi_gridnD(OpPsi_g,Psi_g,Op)
    END IF
  ELSE IF (Basis_IS_allocated(Op%Basis)) THEN
   OpPsi_g(:) = Op%Scalar_g(:) * Psi_g(:)
   OpPsi_g(:) = OpPsi_g(:) -HALF/mass * matmul(Op%Basis%d2gg(:,:,1,1),Psi_g(:))
  END IF

  IF (debug) THEN
    CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
    write(out_unitp,*) 'END OpPsi_grid'
    flush(out_unitp)
  END IF

 END SUBROUTINE OpPsi_grid

 SUBROUTINE calc_OpPsi(Op,Psi,OpPsi)
 USE psi_m, ONLY : psi_t
  TYPE(Op_t),  intent(in)     :: Op
  TYPE(psi_t), intent(in)     :: Psi
  TYPE(psi_t), intent(inout)  :: OpPsi

  IF (.NOT. allocated(Op%RMat)) THEN
    STOP 'ERROR in calc_OpPsi: Op is not initialized!'
  END IF

  IF (allocated(Psi%RVec)) THEN
    OpPsi%RVec = matmul(Op%RMat,Psi%RVec)
  ELSE IF (allocated(Psi%CVec)) THEN
    OpPsi%CVec = matmul(Op%RMat,Psi%CVec)
  ELSE
    STOP 'ERROR in calc_OpPsi: Psi is not initialized!'
  END IF

 END SUBROUTINE calc_OpPsi
END MODULE Op_m
