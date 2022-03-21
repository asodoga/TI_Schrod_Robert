module Op_m
  USE NumParameters_m
  USE diago_m
  USE Basis_m, only : Basis_t
  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis => Null()
    real (kind=Rk),    allocatable :: Scalar_g(:) ! The scalar part of the operator (The potential for Hamiltonian)
    real (kind=Rk),    allocatable :: RMat(:,:)
  END TYPE Op_t

  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi

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

  !  integer :: i

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

    TYPE(Op_t),     intent(inout)       :: Op
    TYPE (Basis_t), intent(in),  target :: Basis
    logical,          parameter         :: debug = .true.
    !logical,         parameter         ::debug = .false.
    integer                             :: ib,iq,iq1,iq2,jb,ib1,ib2
    real (kind=Rk), allocatable         :: Q(:),WT(:),Psi_b(:),Psi_g(:)
    real (kind=Rk), allocatable         :: V(:),OpPsi_g(:),EigenVal(:),EigenVec(:,:)



    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Set_Op'
      flush(out_unitp)
    END IF

    allocate(Psi_b(Basis%nb))
    allocate(Psi_g(Basis%nq))
    allocate(OpPsi_g(Basis%nq))


    CALL alloc_Op(Op,Basis%nb)

    Op%Basis => Basis

    Op%RMat = ZERO

    DO ib=1,Basis%nb
      Psi_b(ib) =ONE
    END DO

    CALL BasisTOGrid_Basis(Psi_g, Psi_b,Basis)

    CALL OpPsi_grid(OpPsi_g,Psi_g,Op)

    ! CALL GridTOBasis_Basis(OpPsi_b, OpPsi_g,Op%Basis)

    DO ib=1,Basis%nb
      Call GridTOBasis_Basis(Op%RMat(:,ib),OpPsi_g,Basis)
    END DO
    CALL write_Op(Op)
    CALL alloc_Op(Op,Basis%nb)

    allocate(EigenVal(Basis%nb))
    allocate(EigenVec(Basis%nb,Basis%nb))
    CALL  diagonalization(Op%RMat,EigenVal,EigenVec,Basis%nb)
    Write(out_unitp,*)
    Write(out_unitp,*)
    Write(out_unitp,*) 'eigenvalues = '
    DO ib=1,Basis%nb
        write(out_unitp,*) EigenVal(ib)
    END DO

    Write(out_unitp,*)
    Write(out_unitp,*)

    DO ib=1,Basis%nb
        write(*,*) (EigenVec(ib,iq),iq=1,Basis%nb)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'END Set_Op'
      flush(out_unitp)
    END IF

  END SUBROUTINE Set_Op

  SUBROUTINE Set_Op1(Op,Basis)
  USE Basis_m
  USE Molec_m

    TYPE(Op_t),     intent(inout)       :: Op
    TYPE (Basis_t), intent(in),  target :: Basis
    logical,          parameter         :: debug = .true.
    !logical,         parameter         ::debug = .false.
    integer                             :: ib,iq,iq1,iq2,jb,ib1,ib2
    real (kind=Rk), allocatable         :: Q(:),WT(:),Psi_b(:),Psi_g(:)
    real (kind=Rk), allocatable         :: V(:),OpPsi_g(:),EigenVal(:),EigenVec(:,:)

    Op%Basis => Basis
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Set_Op'
      flush(out_unitp)
    END IF
    IF(allocated(Basis%tab_basis)) THEN
      allocate(Q(size(Basis%tab_basis)))
      allocate(WT(Basis%nq))
      allocate(Psi_b(Basis%nb))
      allocate(Psi_g(Basis%nq))
      allocate(V(Basis%nq))
      allocate(OpPsi_g(Basis%nq))
      allocate(OP%RMat(Basis%nb,Basis%nb))
      iq = 0
      Do iq1=1,Basis%tab_basis(1)%nq
      DO iq2=1,Basis%tab_basis(2)%nq
        iq=iq+1
        Q(1)=Basis%tab_basis(1)%x(iq1)
        Q(2)=Basis%tab_basis(2)%x(iq2)
        V(iq)=Calc_pot(Q)  !Calc_pot(Q(1))+ Calc_pot(Q(2))
        !V(iq) = Calc_pot(Basis%tab_basis(1)%x(iq1))+Calc_pot(Basis%tab_basis(2)%x(iq2))
        Write(10,*) Q(:),V(iq)
      END DO
      END DO
      ! calculation of Op|b_i>
      Op%RMat = ZERO

      ib1=1
      ib2=0
      DO ib=1,Basis%nb

        IF (ib2 == Basis%tab_basis(2)%nb) THEN
           ib1 = ib1 + 1
           ib2 = 1
         ELSE
           ib2 = ib2 + 1
         END IF

         Psi_b(:)=ZERO
         Psi_b(ib)=ONE

         Call BasisTOGrid_Basis(Psi_g,Psi_b,Basis)

         OpPsi_g(:) = V(:)*Psi_g(:)
         iq1=1
         iq2=0
       DO Iq=1,Basis%nq
         IF(iq2 == Basis%tab_basis(2)%nq) THEN
            iq1 = iq1 + 1
            iq2 = 1
          ELSE
            iq2 = iq2 + 1
          END IF
          OpPsi_g(iq) = OpPsi_g(iq) -HALF/mass *(Basis%tab_basis(1)%d2gb(iq1,ib1,1,1) &
                                               * Basis%tab_basis(2)%d0gb(iq2,ib2))
          OpPsi_g(iq) = OpPsi_g(iq) -HALF/mass *(Basis%tab_basis(1)%d0gb(iq1,ib1) &
                                               * Basis%tab_basis(2)%d2gb(iq2,ib2,1,1))
        END DO
        Call GridTOBasis_Basis(Op%RMat(:,ib),OpPsi_g,Basis)
      END DO
      CALL write_Op(Op)

    ELSE IF( Basis_IS_allocated(Basis)) THEN

      CALL alloc_Op(Op,Basis%nb)
      allocate( Q(1))
      ! calculation of a potential on the grid
      allocate(V(Basis%nq))
      DO iq=1,Basis%nq
        Q = Basis%x(iq)
        V(iq) = Calc_pot(Q)
      END DO

      DO ib=1,Basis%nb
        OpPsi_g = V * Basis%d0gb(:,ib) ! potential part
        OpPsi_g = OpPsi_g -HALF/mass * Basis%d2gb(:,ib,1,1) ! -1/2mass d2./dx2 part
        ! OpPsi_g is a vector on the grid. It must be projected on the basis (integration)
        OpPsi_g = OpPsi_g * Basis%w
        DO jb=1,Basis%nb
           Op%RMat(jb,ib) = dot_product(Basis%d0gb(:,jb),OpPsi_g)
        END DO
      END DO
      CALL write_Op(Op)
      !Test robert
      CALL OpPsi_grid(OpPsi_g,Psi_g,Op)

    ELSE
      STOP 'ERROR in Set_Op: the Basis is not initialized'
    END IF
    allocate(EigenVal(Basis%nb))
    allocate(EigenVec(Basis%nb,Basis%nb))
    CALL  diagonalization(Op%RMat,EigenVal,EigenVec,Basis%nb)
    Write(out_unitp,*)
    Write(out_unitp,*)
    Write(out_unitp,*) 'eigenvalues = '
    DO ib=1,Basis%nb
        write(out_unitp,*) EigenVal(ib)
    END DO

    Write(out_unitp,*)
    Write(out_unitp,*)

    DO ib=1,Basis%nb
        write(*,*) (EigenVec(ib,iq),iq=1,Basis%nb)
    END DO

    IF (debug) THEN
      write(out_unitp,*) 'END '
      flush(out_unitp)
    END IF

  END SUBROUTINE Set_Op1

  SUBROUTINE Potential(V,Basis)
  USE Basis_m
  USE Molec_m
    logical,          parameter         :: debug = .true.
   !logical,         parameter          ::debug = .false.
    integer                             :: iq,iq1,iq2
    real (kind=Rk), allocatable         :: Q(:)
    real (kind=Rk), intent(inout)       :: V(:)
    TYPE(Basis_t),  intent(in)          :: Basis

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Potential'
      flush(out_unitp)
    END IF


    IF (.NOT. Basis_IS_allocated(Basis)) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
      write(out_unitp,*) " the basis is not allocated."
      STOP "ERROR BasisTOGrid_Basis: the basis is not allocated."
    END IF

    IF(allocated(Basis%tab_basis)) THEN
      allocate(Q(size(Basis%tab_basis)))
      iq1=1
      iq2=0
      DO iq=1,Basis%nq

        IF (iq2 == Basis%tab_basis(2)%nq) THEN
          iq1 = iq1 + 1
          iq2 = 1
        ELSE
          iq2 = iq2 + 1
        END IF
        Q(1)=Basis%tab_basis(1)%x(iq1)
        Q(2)=Basis%tab_basis(2)%x(iq2)
        V(iq)=Calc_pot(Q)
     END DO

   ELSE
     allocate( Q(1))
     DO iq=1,Basis%nq
       Q = Basis%x(iq)
       V(iq) = Calc_pot(Q)
     END DO
   END IF

   IF (debug) THEN
    write(out_unitp,*) 'END Potential'
    flush(out_unitp)
   END IF
END SUBROUTINE Potential

SUBROUTINE OpPsi_grid(OpPsi_g,Psi_g,Op)
USE Basis_m
USE Molec_m

    TYPE(Op_t) , intent(inout)          :: Op
    real (kind=Rk), intent(in)          :: Psi_g(:)
    real (kind=Rk), intent(inout)       :: OpPsi_g(:)
    logical,          parameter         :: debug = .true.
   !logical,         parameter         ::debug = .false.
    integer                             :: ib,iq,iq1,iq2,ib1,ib2


    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING OpPsi_grid'
      call Write_op(op)
      flush(out_unitp)
    END IF


    allocate(Op%Scalar_g(Op%Basis%nq))
  stop 'OK'
    CALL Potential(Op%Scalar_g,Op%Basis)

    OpPsi_g(:) = Op%Scalar_g(:)*Psi_g(:)

    OpPsi_g = OpPsi_g -HALF/mass *matmul(Op%Basis%d2gg(:,:,1,1),Psi_g)

    IF (debug) THEN
      write(out_unitp,*) 'END OpPsi_grid'
      flush(out_unitp)
    END IF
  END SUBROUTINE OpPsi_grid

  SUBROUTINE set_op_vo(Op,Basis)
  USE Basis_m
  USE Molec_m
    TYPE(Op_t),     intent(inout)       :: Op
    TYPE (Basis_t), intent(in),  target :: Basis


    integer :: ib,iq,iq1,iq2,jb,ib1,ib2,jb1,jb2
    real (kind=Rk), allocatable :: Q(:),WT(:),G(:),B1(:)
    real (kind=Rk), allocatable :: V(:),OpPsi_g(:),EigenVal(:),EigenVec(:,:)
    IF(allocated(Basis%tab_basis)) THEN
      allocate(Q(size(Basis%tab_basis)))
      allocate(WT(Basis%nq))
      allocate(G(Basis%nq))
      allocate(B1(Basis%nb))
      allocate(V(Basis%nq))
      allocate(OpPsi_g(Basis%nq))
      allocate(OP%RMat(Basis%nb,Basis%nb))
      iq = 0
      Do iq1=1,Basis%tab_basis(1)%nq
      DO iq2=1,Basis%tab_basis(2)%nq
        iq=iq+1
        Q(1)=Basis%tab_basis(1)%x(iq1)
        Q(2)=Basis%tab_basis(2)%x(iq2)
        V(iq)=Calc_pot(Q)  !Calc_pot(Q(1))+ Calc_pot(Q(2))

        !V(iq) = Calc_pot(Basis%tab_basis(1)%x(iq1))+Calc_pot(Basis%tab_basis(2)%x(iq2))
        Write(10,*) Q(:),V(iq)
      END DO
      END DO

      ! calculation of Op|b_i>
      ib=0
      Op%RMat = ZERO
      DO ib1=1,Basis%tab_basis(1)%nb
      DO ib2=1,Basis%tab_basis(2)%nb
         ib=ib+1
         iq=0
         Do iq1=1,Basis%tab_basis(1)%nq
         DO iq2=1,Basis%tab_basis(2)%nq
           iq=iq+1
           OpPsi_g(iq) = V(iq) * Basis%tab_basis(1)%d0gb(iq1,ib1)*Basis%tab_basis(2)%d0gb(iq2,ib2) ! potential part
           OpPsi_g(iq) = OpPsi_g(iq) -HALF/mass *(Basis%tab_basis(1)%d2gb(iq1,ib1,1,1) &
                                                * Basis%tab_basis(2)%d0gb(iq2,ib2))
           OpPsi_g(iq) = OpPsi_g(iq) -HALF/mass *(Basis%tab_basis(1)%d0gb(iq1,ib1) &
                                                * Basis%tab_basis(2)%d2gb(iq2,ib2,1,1))
           WT(iq)=Basis%tab_basis(1)%w(iq1)*Basis%tab_basis(2)%w(iq2)
           OpPsi_g(iq) = OpPsi_g(iq) *WT(iq)
         END DO
         END DO
         jb=0
         DO jb1=1,Basis%tab_basis(1)%nb
         DO jb2=1,Basis%tab_basis(2)%nb
           jb=jb+1
           iq=0
           Do iq1=1,Basis%tab_basis(1)%nq
           DO iq2=1,Basis%tab_basis(2)%nq
             iq=iq+1
             Op%RMat(jb,ib) = Op%RMat(jb,ib)+ Basis%tab_basis(1)%d0gb(iq1,jb1)*Basis%tab_basis(2)%d0gb(iq2,jb2)*OpPsi_g(iq)
           END DO
           END DO
         END DO
         END DO
       END DO
       END DO
       CALL write_Op(Op)
    ELSE IF( Basis_IS_allocated(Basis)) THEN
       CALL alloc_Op(Op,Basis%nb)

      Op%Basis => Basis
      allocate( Q(1))
      allocate(G(Basis%nq))
      allocate(B1(Basis%nb))
      ! calculation of a potential on the grid
      allocate(V(Basis%nq))
      DO iq=1,Basis%nq
        Q = Basis%x(iq)
        V(iq) = Calc_pot(Q)
      END DO
      DO ib=1,Basis%nb
        OpPsi_g = V * Basis%d0gb(:,ib) ! potential part
        OpPsi_g = OpPsi_g -HALF/mass * Basis%d2gb(:,ib,1,1) ! -1/2mass d2./dx2 part
        ! OpPsi_g is a vector on the grid. It must be projected on the basis (integration)
        OpPsi_g = OpPsi_g * Basis%w
        DO jb=1,Basis%nb
           Op%RMat(jb,ib) = dot_product(Basis%d0gb(:,jb),OpPsi_g)
        END DO
      END DO
      CALL write_Op(Op)
    ELSE
      STOP 'ERROR in Set_Op: the Basis is not initialized'
    END IF

    allocate(EigenVal(Basis%nb))
    allocate(EigenVec(Basis%nb,Basis%nb))
    CALL  diagonalization(Op%RMat,EigenVal,EigenVec,Basis%nb)
    Write(out_unitp,*)
    Write(out_unitp,*)
    Write(out_unitp,*) 'eigenvalues = '
    DO ib=1,Basis%nb
        write(out_unitp,*) EigenVal(ib)
    END DO
    Write(out_unitp,*)
    Write(out_unitp,*)

    DO ib=1,Basis%nb
        write(*,*) (EigenVec(ib,iq),iq=1,Basis%nb)
    END DO



  END SUBROUTINE Set_Op_vo

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
