module Op_m
  USE NumParameters_m
  USE diago_m
  USE NDindex_m
  USE Basis_m, only : Basis_t

  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis => Null()
    real (kind=Rk),    allocatable :: Scalar_g(:) ! The scalar part of the operator (The  for Hamiltonian)
    real (kind=Rk),    allocatable :: RMat(:,:)

  END TYPE Op_t

  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi,TEST_OpPsi_grid,Diago_Op,Make_Mat_Op,OpPsi_grid2D

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
    logical,         parameter         ::debug = .false.
    integer                             :: ib,iq
    real (kind=Rk), allocatable         :: Psi_b(:),Psi_g(:)
    real (kind=Rk), allocatable         :: OpPsi_g(:),EigenVal(:),EigenVec(:,:)

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
    real (kind=Rk), allocatable         :: Psi_b(:),Psi_g(:)
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
    DO ib=1,Op%Basis%nb
        write(out_unitp,*) EigenVal(ib)
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
    integer                             :: iq!,iq1,iq2
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
        Q(1)=Op%Basis%tab_basis(1)%x(tab_iq(1))
        Q(2)=Op%Basis%tab_basis(2)%x(tab_iq(2))
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

SUBROUTINE OpPsi_grid3D(OpPsi_g,Psi_g,Op)
USE Basis_m
USE Molec_m
USE UtilLib_m
  TYPE(Op_t) , intent(inout)          :: Op
  real (kind=Rk), intent(in)          :: Psi_g(:)
  real (kind=Rk), intent(inout)       :: OpPsi_g(:)
  real (kind=Rk), allocatable         :: Psi_ggg(:,:,:)
  real (kind=Rk), allocatable         :: opPsi_ggg(:,:,:)
  !logical,          parameter         :: debug = .true.
  logical,         parameter          :: debug = .false.
  logical                             :: Endloop_q
  logical                             :: Endloop_b
  integer                             :: iq,jq1,jq2
  integer,        allocatable         :: tab_iq(:)
  integer,        allocatable         :: tab_iq1(:)
  integer,        allocatable         :: tab_iq2(:)
  integer,        allocatable         :: tab_ib(:)

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING OpPsi_grid3D'
    CALL Write_op(op)
    CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
   flush(out_unitp)
  END IF
  OpPsi_g(:) = Op%Scalar_g(:)*Psi_g(:)

  Allocate(Psi_ggg(Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq,Op%Basis%tab_basis(3)%nq))
  Allocate(OpPsi_ggg(Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq,Op%Basis%tab_basis(3)%nq))

  Psi_ggg(:,:,:) = reshape(Psi_g,shape = [Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq,Op%Basis%tab_basis(3)%nq])
  Oppsi_ggg(:,:,:) = reshape(Oppsi_g,shape = [Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq,Op%Basis%tab_basis(3)%nq])

  Allocate(Tab_iq(size(Op%Basis%tab_basis)))

  CALL Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
  Iq=0
  DO
   Iq=Iq+1
   CALL increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)

   IF (Endloop_q) exit

    OpPsi_ggg(tab_iq(1),tab_iq(2),tab_iq(3)) = OpPsi_ggg(tab_iq(1),tab_iq(2),tab_iq(3))&
     -HALF/mass*dot_product(Op%Basis%tab_basis(1)%d2gg(tab_iq(1),:,1,1),Psi_ggg(tab_iq(1),  :      ,  tab_iq(3)))&
     -HALF/mass*dot_product(Op%Basis%tab_basis(2)%d2gg(tab_iq(2),:,1,1),Psi_ggg(tab_iq(1),tab_iq(2),     :    )) &
     -HALF/mass*dot_product(Op%Basis%tab_basis(3)%d2gg(tab_iq(3),:,1,1),Psi_ggg(   :     ,tab_iq(2),  tab_iq(3)))
  END DO
  Deallocate(Tab_iq)

  Oppsi_g(:) = reshape(Oppsi_ggg,shape=[Op%Basis%nq])
  
  IF (debug) THEN
  	CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
  	write(out_unitp,*) 'END OpPsi_grid3D'
  	flush(out_unitp)
  END IF

END SUBROUTINE OpPsi_grid3D

SUBROUTINE OpPsi_grid2D(OpPsi_g,Psi_g,Op)
USE Basis_m
USE Molec_m
USE UtilLib_m
    TYPE(Op_t) , intent(inout)          :: Op
    real (kind=Rk), intent(in)          :: Psi_g(:)
    real (kind=Rk), intent(inout)       :: OpPsi_g(:)
    real (kind=Rk), allocatable         :: Psi_gg(:,:)
    real (kind=Rk), allocatable         :: opPsi_gg(:,:)
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.
    logical                             :: Endloop_q
    logical                             :: Endloop_b
    integer                             :: iq,jq1,jq2
    integer,        allocatable         :: tab_iq(:)
    integer,        allocatable         :: tab_iq1(:)
    integer,        allocatable         :: tab_ib(:)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING OpPsi_grid2d'
      CALL Write_op(op)
      CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
      flush(out_unitp)
    END IF
    OpPsi_g(:) = Op%Scalar_g(:)*Psi_g(:)
    Allocate(Psi_gg(Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq))
    Allocate(OpPsi_gg(Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq))

    Psi_gg(:,:) = reshape(Psi_g,shape = [Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq])
    Oppsi_gg(:,:) = reshape(Oppsi_g,shape = [Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq])

    Allocate(Tab_iq(size(Op%Basis%tab_basis)))

    Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
    Iq=0
    Do
      Iq=Iq+1
      CALL increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)

      IF (Endloop_q) exit

      OpPsi_gg(tab_iq(1),tab_iq(2)) = OpPsi_gg(tab_iq(1),tab_iq(2))&
      -HALF/mass*dot_product(Op%Basis%tab_basis(1)%d2gg(tab_iq(1),:,1,1),Psi_gg(:,tab_iq(2)))&
      -HALF/mass*dot_product(Op%Basis%tab_basis(2)%d2gg(tab_iq(2),:,1,1),Psi_gg(tab_iq(1),:))
    END DO
    Deallocate(Tab_iq)

    Oppsi_g(:) = reshape(Oppsi_gg,shape=[Op%Basis%nq])
    IF (debug) THEN
    	CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
    	write(out_unitp,*) 'END OpPsi_grid2D'
    	flush(out_unitp)
    END IF

END SUBROUTINE OpPsi_grid2D

SUBROUTINE OpPsi_grid2D0(OpPsi_g,Psi_g,Op)
USE Basis_m
USE Molec_m
USE UtilLib_m
    TYPE(Op_t) , intent(inout)          :: Op
    real (kind=Rk), intent(in)          :: Psi_g(:)
    real (kind=Rk), intent(inout)       :: OpPsi_g(:)
    real (kind=Rk), allocatable         :: Psi_gg(:,:)
    real (kind=Rk), allocatable         :: opPsi_gg(:,:)
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.
    logical                             :: Endloop_q
    logical                             :: Endloop_b
    integer                             :: iq,jq1,jq2
    integer,        allocatable         :: tab_iq(:)
    integer,        allocatable         :: tab_iq1(:)
    integer,        allocatable         :: tab_ib(:)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING OpPsi_grid2d'
      CALL Write_op(op)
      CALL Write_RVec(Psi_g,out_unitp,5,name_info='Psi_g')
      flush(out_unitp)
    END IF
      OpPsi_g(:) = Op%Scalar_g(:)*Psi_g(:)
      Allocate(Psi_gg(Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq))
      Allocate(OpPsi_gg(Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq))

      Psi_gg(:,:) = reshape(Psi_g,shape = [Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq])
      Oppsi_gg(:,:) = reshape(Oppsi_g,shape = [Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq])

      Allocate(Tab_iq(size(Op%Basis%tab_basis)))

      Call Init_tab_ind(Tab_iq,Op%Basis%NDindexq)
      Iq=0
      Do
        Iq=Iq+1
        CALL increase_NDindex(Tab_iq,Op%Basis%NDindexq,Endloop_q)

        IF (Endloop_q) exit

        OpPsi_gg(tab_iq(1),tab_iq(2)) = OpPsi_gg(tab_iq(1),tab_iq(2))-HALF/mass &
        *dot_product(Op%Basis%tab_basis(1)%d2gg(tab_iq(1),:,1,1),Psi_gg(:,tab_iq(2)))

      END DO
      Deallocate(Tab_iq)
      Allocate(Tab_iq1(size(Op%Basis%tab_basis)))
      Call Init_tab_ind(Tab_iq1,Op%Basis%NDindexq)
      Iq=0
      DO
        Iq=Iq+1
        CALL increase_NDindex(Tab_iq1,Op%Basis%NDindexq,Endloop_q)
        IF (Endloop_q) exit

        OpPsi_gg(tab_iq1(1),tab_iq1(2)) =OpPsi_gg(tab_iq1(1),tab_iq1(2))-HALF/mass &
        *dot_product(Op%Basis%tab_basis(2)%d2gg(tab_iq1(2),:,1,1),Psi_gg(tab_iq1(1),:))

      END DO

      Deallocate(Tab_iq1)
      Oppsi_g(:) = reshape(Oppsi_gg,shape=[Op%Basis%nq])
      IF (debug) THEN
      	CALL Write_RVec(OpPsi_g,out_unitp,5,name_info='OpPsi_g')
      	write(out_unitp,*) 'END OpPsi_grid2D'
      	flush(out_unitp)
      END IF

END SUBROUTINE OpPsi_grid2D0
SUBROUTINE OpPsi_grid(OpPsi_g,Psi_g,Op)
USE Basis_m
USE Molec_m
USE UtilLib_m

    TYPE(Op_t) , intent(inout)          :: Op
    real (kind=Rk), intent(in)          :: Psi_g(:)
    real (kind=Rk), intent(inout)       :: OpPsi_g(:)
    !real (kind=Rk), allocatable         :: Psi_gg(:,:)
    !real (kind=Rk), allocatable         :: opPsi_gg(:,:)
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




      IF(Size(Op%Basis%tab_basis)==3) THEN

        CALL OpPsi_grid3D(OpPsi_g,Psi_g,Op)

      ELSE IF(Size(Op%Basis%tab_basis)==2)  THEN

        CALL OpPsi_grid2D(OpPsi_g,Psi_g,Op)

      ELSE

        STOP ' OpPsi_grid STOPED for Overflow  '

      END IF

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





  SUBROUTINE TEST_OpPsi_grid(Op)
  USE Basis_m
  USE Molec_m
    TYPE(Op_t) , intent(inout)          :: Op
    real (kind=Rk),allocatable          :: Psi_g(:)
    real (kind=Rk),allocatable          :: Psi_gg(:,:)
    real (kind=Rk),allocatable          :: Dif_g(:)
    real (kind=Rk),allocatable          :: OpPsi_g(:)
    !logical,          parameter         :: debug = .true.
    logical,         parameter          :: debug = .false.
    integer                             :: iq,iq1,iq2
    real(kind=Rk),    parameter         :: eps = TEN**(-TEN)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING TEST_CALC'
      CALL write_op(OP)
      flush(out_unitp)
    END IF

    IF(allocated(Op%Basis%tab_basis)) THEN

      allocate( Psi_g(Op%Basis%nq))
      allocate( Psi_gg(Op%Basis%tab_basis(1)%nq,Op%Basis%tab_basis(2)%nq))
      ! Psi analytic
        DO iq1=1,Op%Basis%tab_basis(1)%nq
        DO iq2=1,Op%Basis%tab_basis(2)%nq

           Psi_gg(iq1,iq2)=exp(-HALF*(Op%Basis%tab_basis(1)%x(iq1))**2)*exp(-HALF*(Op%Basis%tab_basis(2)%x(iq2))**2)

        END DO
        END DO
       Psi_g(:) = reshape(Psi_gg,shape=[Op%Basis%nq])
       allocate( OpPsi_g(Op%Basis%nq))

      CALL OpPsi_grid(OpPsi_g,Psi_g,Op)

      DO iq=1,Op%Basis%nq
        Write(out_unitp,*)iq,OpPsi_g(iq),ONE*Psi_g(iq)
      END DO
      allocate( Dif_g(Op%Basis%nq))

      DO iq=1,Op%Basis%nq
        Dif_g(iq)=ABS(OpPsi_g(iq)-ONE*Psi_g(iq))
      END DO

      IF (eps.gt.maxval(Dif_g(:)) ) THEN
        Write(out_unitp,*) 'OpPsi_grid is correct '
        Write(out_unitp,*) maxval(Dif_g(:))
      ELSE
        Write(out_unitp,*) maxval(Dif_g(:))
        STOP 'OpPsi_grid is not correct'
      END IF

    ELSE IF( Basis_IS_allocated(Op%Basis)) THEN

      allocate( Psi_g(Op%Basis%nq))

     ! Psi analytic

      DO iq=1,Op%Basis%nq
        Psi_g(iq)=exp(-HALF*(Op%Basis%x(iq))**TWO)
      END DO

      allocate( OpPsi_g(Op%Basis%nq))

      CALL OpPsi_grid(OpPsi_g,Psi_g,Op)
      DO iq=1,Op%Basis%nq
        Write(out_unitp,*)iq,OpPsi_g(iq),HALF*Psi_g(iq)
      END DO
      allocate( Dif_g(Op%Basis%nq))

      DO iq=1,Op%Basis%nq
        Dif_g(iq)=ABS(OpPsi_g(iq)-HALF*Psi_g(iq))
      END DO

      IF (eps.gt.maxval(Dif_g(:)) ) THEN
        Write(out_unitp,*) 'OpPsi_grid is correct '
        Write(out_unitp,*) maxval(Dif_g(:))
      ELSE
        Write(out_unitp,*) maxval(Dif_g(:))
        STOP 'OpPsi_grid is not correct'
      END IF

    END IF

    IF (debug) THEN
    !  CALL Write_RVec(Dif_g,out_unitp,5,name_info='Dif_g')
      write(out_unitp,*) 'END TEST_CALC'
      flush(out_unitp)
    END IF

  END SUBROUTINE TEST_OpPsi_grid

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
