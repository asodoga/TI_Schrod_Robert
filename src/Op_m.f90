module Op_m
  USE NumParameters_m
  USE diago_m
  USE NDindex_m
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


  SUBROUTINE Set_Op(Op,Basis)
  USE Basis_m
  USE Molec_m
  USE NDindex_m

    TYPE(Op_t),     intent(inout)       :: Op
    TYPE (Basis_t), intent(in),  target :: Basis
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.
    integer                             :: ib,iq


    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Set_Op'
      call write_basis(Basis)
      flush(out_unitp)
    END IF


    Op%Basis => Basis

    IF (debug) THEN
      write(out_unitp,*) 'END Set_Op'
      flush(out_unitp)
    END IF

  END SUBROUTINE Set_Op

  SUBROUTINE Make_Mat_OP(Op)
  USE Basis_m
    TYPE (Op_t),     intent(inout)      :: Op
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.
    integer                             :: ib,iq,jb
    real (kind=Rk), pointer             :: Psi_g(:) => null()
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
    deallocate(Psi_b)
    deallocate(Psi_g)
    deallocate(OpPsi_g)

    IF (debug) THEN
      write(out_unitp,*) 'END Make_Mat_OP'
      flush(out_unitp)
    END IF

  END SUBROUTINE Make_Mat_OP

  SUBROUTINE Diago_Op(Op)
  USE Basis_m
  USE Molec_m

    TYPE(Op_t),     intent(inout)       :: Op
    !logical,          parameter         :: debug = .true.
    logical,         parameter          ::debug = .false.
    integer                             :: ib,jb
    real (kind=Rk), allocatable         :: EigenVal(:),EigenVec(:,:)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Diago_Op'
      flush(out_unitp)
    END IF

    allocate(EigenVal(Op%Basis%nb))
    allocate(EigenVec(Op%Basis%nb,Op%Basis%nb))
    CALL  diagonalization(Op%RMat,EigenVal,EigenVec,Op%Basis%nb)
    Write(out_unitp,*)
    Write(out_unitp,*)
    Write(out_unitp,*) 'eigenvalues = '
    DO ib=1,3!Op%Basis%nb
        write(out_unitp,*) EigenVal(ib)*8066_RK*27.21_RK
    END DO
    Write(out_unitp,*)
    Write(out_unitp,*)
    DO ib=1,Op%Basis%nb
        write(*,*) (EigenVec(ib,jb),jb=1,Op%Basis%nb)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'END Diago_Op'
      flush(out_unitp)
    END IF
  END SUBROUTINE Diago_Op

  SUBROUTINE Potential(Op)
   USE Basis_m
   USE Molec_m
   !logical,          parameter         :: debug = .true.
   logical,         parameter          :: debug = .false.
   logical                             :: Endloop_q
   integer ,allocatable                :: tab_iq(:)
   integer                             :: iq,inb
   real (kind=Rk), allocatable         :: Q(:)
   TYPE(Op_t),  intent(inout)          :: Op

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
     Iq=0
     DO
      Iq=Iq+1
      CALL increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)
      IF (Endloop_q) exit
      DO inb=1, size(Op%Basis%tab_basis)
        Q(inb)=Op%Basis%tab_basis(inb)%x(tab_iq(inb))
      END DO
      Op%Scalar_g(iq)=Calc_pot(Q)
     END DO
     Deallocate(Tab_iq)
     Deallocate(Q)
   ELSE
     allocate( Q(1))
     DO iq=1,Op%Basis%nq
       Q = Op%Basis%x(iq)
       Op%Scalar_g(iq) = Calc_pot(Q)
     END DO
   END IF

   IF (debug) THEN
     write(out_unitp,*)'v', Op%Scalar_g
     write(out_unitp,*) 'END Potential'
     flush(out_unitp)
   END IF
 END SUBROUTINE Potential





 !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ Robert en problème §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!!!!!!!!!
  SUBROUTINE KEOiPsi_grid(OpPsi_g,Psi_g,Op,Op%Grid(iterm)%Vec,inb)
  USE Basis_m
  USE Molec_m
  USE UtilLib_m
    TYPE(Op_t),     intent(in),target   :: Op
    real (kind=Rk), intent(in) ,target  :: Psi_g(:)
    real (kind=Rk), intent(inout),target:: OpPsi_g(:)
    integer,        intent(in)          :: inb
    real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
    real (kind=Rk), pointer             :: d2gg(:,:)
    real (kind=Rk), pointer             :: d1gg(:,:)
    real (kind=Rk), pointer             :: KEO1psi_g(:)
    real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.
    integer                             :: iq,i1,i3,inb
    integer                             :: Iq1,Iq2,Iq3

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING KEOiPsi_grid'
      CALL Write_op(op)
      CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
     flush(out_unitp)
    END IF
    inb = Op%Grid(2)%DerivIndex(1)
    !inb2 = Op%Grid(2)%DerivIndex(2) = 0

    !!!!1) Action de d/dQi : k1_psi = d/dQi psi

    Iq1 = Product(Op%Basis%tab_basis(1:inb-1)%nq)

    Iq2 = Op%Basis%tab_basis(inb)%nq

    Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inb+1:-1)%nq)

    OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO1Psi_g

    Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g

    d1gg(1:Op%Basis%tab_basis(inb)%nq,1:Op%Basis%tab_basis(inb)%nq)=>Op%Basis%tab_basis(inb)%d1gg

    DO i3=1,ubound(Psi_ggg,dim=3)
    DO i1=1,ubound(Psi_ggg,dim=1)

      OpPsi_ggg(i1,:,i3) = OpPsi_ggg(i1,:,i3)+ matmul(d1gg,Psi_ggg(i1,:,i3))

    END DO
    END DO
    !3) action de f1^j(Q) :  kpsi = f1^j(Q) * k1_psi
    OpPsi_g(:) =   Op%Grid(2)%Vec(:)* KEO1Psi_g(:)

    IF (debug) THEN
    	CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
    	write(out_unitp,*) 'END KEOiPsi_grid'
    	flush(out_unitp)
    END IF

  END SUBROUTINE KEOiPsi_grid


!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§ Robert en problème §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!!!!!!!!!

 SUBROUTINE KEOijPsi_grid(OpPsi_g,Psi_g,Op,Op%Grid(iterm)%Vec,inb1,inb2)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   real (kind=Rk), intent(in) ,target  :: Psi_g(:)
   real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
   real (kind=Rk), pointer             :: d2gg(:,:)
   real (kind=Rk), pointer             :: d1gg(:,:)
   real (kind=Rk), pointer             :: KEO1psi_g(:)
   real (kind=Rk), pointer             :: KEO2psi_g(:)
   real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
   !logical,          parameter         :: debug = .true.
   logical,         parameter          :: debug = .false.
   integer                             :: iq,i1,i3,inb
   integer                             :: Iq1,Iq2,Iq3

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
     CALL Write_op(op)
     CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
   END IF

   !!!!1) Action de d/dQi : k1_psi = d/dQi psi

   Iq1 = Product(Op%Basis%tab_basis(1:inb1-1)%nq)

   Iq2 = Op%Basis%tab_basis(inb1)%nq

   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inb1+1:-1)%nq)

   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO1Psi_g

   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g

   d1gg(1:Op%Basis%tab_basis(inb1)%nq,1:Op%Basis%tab_basis(inb1)%nq)=>Op%Basis%tab_basis(inb1)%d1gg

   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)

     OpPsi_ggg(i1,:,i3) = OpPsi_ggg(i1,:,i3)+ matmul(d1gg,Psi_ggg(i1,:,i3))

   END DO
   END DO

 !!!!2) Action de d/dQj :   k2_psi = d/dQj k1_psi


   Iq1 = Product(Op%Basis%tab_basis(1:inb2-1)%nq)

   Iq2 = Op%Basis%tab_basis(inb2)%nq

   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inb2+1:-1)%nq)

   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO2Psi_g

   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => KEO1Psi_g

   d1gg(1:Op%Basis%tab_basis(inb2)%nq,1:Op%Basis%tab_basis(inb2)%nq)=>Op%Basis%tab_basis(inb2)%d1gg

   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)

     OpPsi_ggg(i1,:,i3) = OpPsi_ggg(i1,:,i3) + matmul(d1gg,Psi_ggg(i1,:,i3))

   END DO
   END DO


  !3) action de f2^ij(Q) :  kpsi = f2^ij(Q) * k2_psi

   OpPsi_g(:) =   Op%Grid(iterm)%Vec(:)*KEO2Psi_g(:)


   IF (debug) THEN
   	CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
   	write(out_unitp,*) 'END KEOijPsi_grid'
   	flush(out_unitp)
   END IF

 END SUBROUTINE KEOijPsi_grid

 SUBROUTINE KEOiiPsi_grid(OpPsi_g,Psi_g,Op,Op%Grid(iterm)%Vec,inb)
 USE Basis_m
 USE Molec_m
 USE UtilLib_m
   TYPE(Op_t),     intent(in),target   :: Op
   real (kind=Rk), intent(in) ,target  :: Psi_g(:)
   real (kind=Rk), intent(inout),target:: OpPsi_g(:)
   real (kind=Rk), pointer             :: Psi_ggg(:,:,:)
   real (kind=Rk), pointer             :: d2gg(:,:)
   real (kind=Rk), pointer             :: d1gg(:,:)
   real (kind=Rk), pointer             :: KEO1psi_g(:)
   real (kind=Rk), pointer             :: KEO2psi_g(:)
   real (kind=Rk), pointer             :: OpPsi_ggg(:,:,:)
   !logical,          parameter         :: debug = .true.
   logical,         parameter          :: debug = .false.
   integer                             :: iq,i1,i3,inb
   integer                             :: Iq1,Iq2,Iq3

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING KEOijPsi_grid'
     CALL Write_op(op)
     CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
    flush(out_unitp)
   END IF

   !!!!1) Action de d/dQi : k1_psi = d/dQi psi

   Iq1 = Product(Op%Basis%tab_basis(1:inb-1)%nq)

   Iq2 = Op%Basis%tab_basis(inb)%nq

   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inb+1:-1)%nq)

   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO1Psi_g

   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => Psi_g

   d1gg(1:Op%Basis%tab_basis(inb)%nq,1:Op%Basis%tab_basis(inb)%nq)=>Op%Basis%tab_basis(inb)%d1gg

   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)

     OpPsi_ggg(i1,:,i3) = OpPsi_ggg(i1,:,i3)+ matmul(d1gg,Psi_ggg(i1,:,i3))

   END DO
   END DO

 !!!!2) Action de d/dQj :   k2_psi = d/dQj k1_psi


   Iq1 = Product(Op%Basis%tab_basis(1:inb-1)%nq)

   Iq2 = Op%Basis%tab_basis(inb)%nq

   Iq3 = Product(Op%Basis%tab_basis(size(Op%Basis%tab_basis):inb+1:-1)%nq)

   OpPsi_ggg(1:Iq1,1:Iq2,1:iq3)=> KEO2Psi_g

   Psi_ggg(1:Iq1,1:Iq2,1:iq3)  => KEO1Psi_g

   d1gg(1:Op%Basis%tab_basis(inb)%nq,1:Op%Basis%tab_basis(inb)%nq)=>Op%Basis%tab_basis(inb)%d1gg

   DO i3=1,ubound(Psi_ggg,dim=3)
   DO i1=1,ubound(Psi_ggg,dim=1)

     OpPsi_ggg(i1,:,i3) = OpPsi_ggg(i1,:,i3) + matmul(d1gg,Psi_ggg(i1,:,i3))

   END DO
   END DO


  !3) action de f2^ij(Q) :  kpsi = f2^ij(Q) * k2_psi

   OpPsi_g(:) =   Op%Grid(iterm)%Vec(:)*KEO2Psi_g(:)


   IF (debug) THEN
     CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     write(out_unitp,*) 'END KEOiiPsi_grid'
     flush(out_unitp)
   END IF

 END SUBROUTINE KEOiiPsi_grid



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
    integer                             :: iq,i1,i3,inb1,inb2
    integer                             :: Iq1,Iq2,Iq3


     IF (debug) THEN
       write(out_unitp,*) 'BEGINNING OpPsi_gridcnD'
       CALL Write_op(op)
       CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
      flush(out_unitp)
     END IF


     OpPsi_g(:) = Op%Scalar_g(:)*Psi_g(:)


     !DO inb = 1,size(Op%Basis%tab_basis)
     ! 2e façon
     Opsi(:) = 0
     DO iterm=1,size(Op%Grid)

       inb1 = Op%Grid(iterm)%DerivIndex(1) ! pas compris
       inb2 = Op%Grid(iterm)%DerivIndex(2)

       IF (inb1 > 0 .and. inb2 > 0 .and. inb1 /= inb2) THEN ! terme croisé

         CALL KEOijPsi_grid(OpPsi_g,Psi_g,Op,Op%Grid(iterm)%Vec,inb1,inb2)

       ELSE IF (inb1 > 0 .and. inb2 > 0 .and. inb1 == inb2) THEN ! terme non croisé

         CALL KEOiiPsi_grid(OpPsi_g,Psi_g,Op,Op%Grid(iterm)%Vec,inb1)

       ELSE IF (inb1 > 0 .and. inb2 == 0) THEN

         CALL KEOiPsi_grid(OpPsi_g,Psi_g,Op,Op%Grid(iterm)%Vec,inb1) ! dérivé premiere

       END IF
     END DO

     !END DO

     IF (debug) THEN
     	CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
     	write(out_unitp,*) 'END OpPsi_gridcnD'
     	flush(out_unitp)
     END IF

   END SUBROUTINE OpPsi_gridcnD
! §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!
!§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!!!!!!!!!

 SUBROUTINE OpPsi_gridnD(OpPsi_g,Psi_g,Op)
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
  ! §§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!

  !§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§!!!!!!!!!
  SUBROUTINE OpPsi_grid(OpPsi_g,Psi_g,Op)
  USE Basis_m
  USE Molec_m
  USE UtilLib_m

    TYPE(Op_t) , intent(inout)          :: Op
    real (kind=Rk), intent(in) ,target  :: Psi_g(:)
    real (kind=Rk), intent(inout)       :: OpPsi_g(:)
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.


    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING OpPsi_grid'
      CALL Write_op(op)
      CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
      flush(out_unitp)
    END IF

    IF (.NOT.allocated(Op%Scalar_g)) THEN

      allocate(Op%Scalar_g(Op%Basis%nq))

      CALL Potential(Op)

    END IF


    IF(allocated(Op%Basis%tab_basis)) THEN

     CALL OpPsi_gridnD(OpPsi_g,Psi_g,Op)

    ELSE IF( Basis_IS_allocated(Op%Basis)) THEN

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
