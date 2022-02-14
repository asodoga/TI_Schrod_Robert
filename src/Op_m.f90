module Op_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis

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

  END SUBROUTINE dealloc_Op

  SUBROUTINE write_Op(Op)
    USE UtilLib_m, only : Write_RMat

    TYPE(Op_t), intent(in) :: Op

    integer :: i

    IF (associated(Op%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to Op.'
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


    integer :: ib,iq,iq1,iq2,jb,ib1,ib2,jb1,jb2
    real (kind=Rk), allocatable :: Q(:)
    real (kind=Rk), allocatable :: V(:),OpPsi_g(:)

    IF(allocated(Basis%tab_basis)) THEN
      allocate(Q(size(Basis%tab_basis)))
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
      !  Write(10,*) "iq1,iq2,iq,V(iq)" ,iq1,iq2,iq,Q(:),V(iq)
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
         END DO
         END DO
    !    OpPsi_g = OpPsi_g -HALF/mass * Basis%tab_basis(1)%d2gb(:,ib,1,1) -HALF/mass * Basis%tab_basis(2)%d2gb(:,ib,1,1)
        ! OpPsi_g is a vector on the grid. It must be projected on the basis (integration)
      !  OpPsi_g = OpPsi_g * Basis%tab_basis(1)%w*Basis%tab_basis(2)%w
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
            STOP
    ELSE IF( Basis_IS_allocated(Basis)) THEN
        CALL alloc_Op(Op,Basis%nb)

      Op%Basis => Basis
      allocate( Q(1))

      ! calculation of a potential on the grid
      allocate(V(Basis%nq))
      DO iq=1,Basis%nq
        Q = Basis%x(iq)
        V(iq) = Calc_pot(Q)
      END DO

      ! calculation of Op|b_i>
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
    !END IF

  !  IF (.NOT. Basis_IS_allocated(Basis)) THEN
      STOP 'ERROR in Set_Op: the Basis is not initialized'
    END IF



  END SUBROUTINE Set_Op


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

end module Op_m
