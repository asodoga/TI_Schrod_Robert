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
MODULE Basis_m
  USE NumParameters_m
  USE UtilLib_m
  USE NDindex_m
  USE UtilMath_m


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Basis_t,Basis_IS_Allocated,BasisTOGrid_Basis,GridTOBasis_Basis
  PUBLIC :: Test_Passage,Calc_dngg_grid,Basis_IS_Allocatedtot,write_basis
  PUBLIC :: BasisTOGrid_Basis_rapide,GridTOBasis_Basis_rapide,Read_Construct_Basis
  PUBLIC :: BasisTOGrid_Basis_rapide1,GridTOBasis_Basis_rapide1,NDSmolyak_t
  !PUBLIC :: Read_Basis_old!,Calc_n_smol


  TYPE :: NDSmolyak_t
    Integer, allocatable         :: DSmol(:)
  END TYPE NDSmolyak_t

  TYPE :: Basis_t
    Integer                      :: nb_basis   = 0
    Integer                      :: nb         = 0
    Integer                      :: nq         = 0
    Integer                      :: ndim       = 1
    Integer                      :: A_smol, B_smol
    Integer, allocatable         :: Tags1(:)
    Integer, allocatable         :: Tab_icord_to_ib(:)
    Integer, allocatable         :: Tab_icord_to_idim(:)
    Integer, allocatable         :: Ind_map(:)
    Integer, allocatable         :: tab_ic_to_ib(:)
    Integer, allocatable         :: tab_ic_to_i_dim(:)
    Integer, allocatable         :: Ind_iq(:)
    Integer, allocatable         :: Nbi_li(:,:)
    Real(kind=Rk)                :: A,B,scaleQ,Q0
    Character(len=:),allocatable :: Basis_name
    Real(kind=Rk),   allocatable :: x(:,:)
    Real(kind=Rk),   allocatable :: w(:)
    Real(kind=Rk),   allocatable :: d0gb(:,:)      ! basis s d0gb(nq,nb)
    Real(kind=Rk),   allocatable :: d1gb(:,:,:)    ! basis s d2gb(nq,nb,1)
    Real(kind=Rk),   allocatable :: d1gg(:,:,:)    ! basis s d2gg(nq,nq,1)
    Real(kind=Rk),   allocatable :: d2gb(:,:,:,:)  ! basis s d2gb(nq,nb,1,1)
    Real(kind=Rk),   allocatable :: d2gg(:,:,:,:)  ! basis s d2gg(nq,nq,1,1)
    TYPE(NDindex_t)              :: NDindexq
    TYPE(NDindex_t)              :: NDindexb
    TYPE(NDindex_t)              :: NDindexl
    TYPE(NDSmolyak_t)            :: Smolyak
    TYPE(Basis_t),  allocatable  :: tab_basis(:)
    TYPE(Basis_t),  allocatable  :: tab_Smolyak(:)
    TYPE(Basis_t),  allocatable  :: tab_in_Smolyak(:)
  END TYPE Basis_t

CONTAINS
  RECURSIVE FUNCTION Basis_IS_Allocated(Basis) RESULT(alloc)
    TYPE(Basis_t),   intent(in)  :: Basis
    Logical                      :: alloc
    Integer                      :: i

    alloc = Allocated(Basis%tab_basis)
    IF ( Allocated(Basis%tab_basis)) THEN
      Do i=1,size(Basis%tab_basis)
        alloc  = alloc .and. Basis_IS_Allocated(Basis%tab_basis(i))
      END DO
    ELSE
      alloc =             Allocated(Basis%x)
      alloc = alloc .AND. Allocated(Basis%w)
      alloc = alloc .AND. Allocated(Basis%d0gb)
      alloc = alloc .AND. Allocated(Basis%d1gb)
      alloc = alloc .AND. Allocated(Basis%d2gb)
    END IF
  END FUNCTION Basis_IS_Allocated

  RECURSIVE FUNCTION Basis_IS_Allocatedtot(Basis) RESULT(alloc)

      TYPE(Basis_t),   intent(in)  :: Basis
      Logical                      :: alloc
      Integer                      :: i

      alloc = Allocated(Basis%tab_basis)
      IF ( Allocated(Basis%tab_basis)) THEN
        Do i=1,size(Basis%tab_basis)
          alloc  = alloc .and. Basis_IS_Allocated(Basis%tab_basis(i))
        END DO
      ELSE
        alloc =             Allocated(Basis%x)
        alloc = alloc .AND. Allocated(Basis%w)
        alloc = alloc .AND. Allocated(Basis%d0gb)
        alloc = alloc .AND. Allocated(Basis%d1gb)
        alloc = alloc .AND. Allocated(Basis%d2gb)
        alloc = alloc .AND. Allocated(Basis%d1gg)
        alloc = alloc .AND. Allocated(Basis%d2gg)
      END IF

  END FUNCTION Basis_IS_Allocatedtot

  RECURSIVE SUBROUTINE Write_Basis(Basis)
   USE UtilLib_m
   TYPE(Basis_t),       intent(in)  :: Basis
   Integer                          :: i

   write(out_unitp,*) '-------------------------------------------------'
   write(out_unitp,*) 'Write_Basis'
   write(out_unitp,*) 'nb,nq',Basis%nb,Basis%nq
  ! If(Basis%ndim==1); size(Basis%tab_basis)=1
   write(out_unitp,*) 'ndim','size_tab',Basis%ndim,size(Basis%tab_basis)

     write(out_unitp,*)"A et B",Basis%A,Basis%B
     write(out_unitp,*) 'Scaleq',Basis%Scaleq
     write(out_unitp,*) 'Q0' ,Basis%Q0
     IF (.NOT.Allocated(Basis%x)) THEN
      write(out_unitp,*)' Basis table x is not Allocated.'
     ELSE
      Do i=1,Basis%ndim
        Call Write_RVec(Basis%x(:,i),out_unitp,5,name_info='x')
      END DO
      ! Call Write_RVec(Basis%x(:,2),out_unitp,5,name_info='x')
     END IF
     write(out_unitp,*)
     IF (.NOT.Allocated(Basis%W)) THEN
       write(out_unitp,*)' Basis table w is not Allocated.'
     ELSE

       Call Write_RVec(Basis%w(:),out_unitp,5,name_info='w')
     END IF
     write(out_unitp,*)
     IF (.NOT.Allocated(Basis%d0gb)) THEN
       write(out_unitp,*)' Basis table d0gb is not Allocated.'
     ELSE
       Call Write_RMat(Basis%d0gb,out_unitp,5,name_info='d0gb')
     END IF
     write(out_unitp,*)
     IF (.NOT.Allocated(Basis%d1gb)) THEN
       write(out_unitp,*)' Basis table d1gb is not Allocated.'
     ELSE
       Call Write_RMat(Basis%d1gb(:,:,1),out_unitp,5,name_info='d1gb')
     END IF
     write(out_unitp,*)
     IF (.NOT.Allocated(Basis%d1gg)) THEN
       write(out_unitp,*)' Basis table d1gb is not Allocated.'
     ELSE
       Call Write_RMat(Basis%d1gg(:,:,1),out_unitp,5,name_info='d1gg')
     END IF
     write(out_unitp,*)
     IF (.NOT.Allocated(Basis%d2gb)) THEN
       write(out_unitp,*)' Basis table d1gb is not Allocated.'
     ELSE
       Call Write_RMat(Basis%d2gb(:,:,1,1),out_unitp,5,name_info='d2gb')
     END IF
     write(out_unitp,*)
     IF (.NOT.Allocated(Basis%d2gg)) THEN
       write(out_unitp,*)' Basis table d2gg is not Allocated.'
     ELSE
       Call Write_RMat(Basis%d2gg(:,:,1,1),out_unitp,5,name_info='d2gg')
     END IF

   !  write(out_unitp,*) 'nb_basis',Basis%nb_basis
   IF (Allocated(Basis%tab_basis)) THEN
     DO i=1,size(Basis%tab_basis)
       Call Write_Basis(Basis%tab_basis(i))
     END DO
   END IF
   write(out_unitp,*) '-------------------------------------------------'

 END SUBROUTINE Write_Basis

  SUBROUTINE Read_Basis_simple(name,nb_basis,LB,A_smol,B_smol,nb,nq,A,B,scaleQ,Q0,nio)
   USE UtilLib_m
    !Logical,        parameter               :: debug = .true.
    Logical,       parameter                :: debug = .false.
    Integer,        intent(in)              :: nio
    Integer ,       intent(inout)           :: nb,nq,nb_basis,LB
    Integer ,       intent(inout)           :: A_smol,B_smol
    Integer                                 :: err_io
    Character (len=Name_len) ,intent(inout) :: name
    Real(kind=Rk) ,intent(inout)            :: A,B,scaleQ,Q0

    NAMELIST /basis_nD/ name,nb_basis,LB,nb,nq,A_smol,B_smol,A,B,scaleQ,Q0

    nb_basis  = 0
    nb        = 0
    nq        = 0
    LB        = 0
    A_smol    = 0
    B_smol    = 0
    A         = ZERO
    B         = ZERO
    Q0        = ZERO
    scaleQ    = ONE
    name      = '0'

    read(nio,nml=basis_nD,IOSTAT=err_io)
    write(out_unitp,nml=basis_nD)
  END SUBROUTINE Read_Basis_simple

  SUBROUTINE Read_Construct_Basis(Basis,nio)
   USE UtilLib_m
   IMPLICIT NONE
    !Logical,             parameter      :: debug = .true.
    Logical,             parameter       :: debug = .false.
    TYPE(Basis_t),       intent(inout)  :: Basis
    Integer,             intent(in)     :: nio
    Integer, allocatable                :: NDend_q(:)
    Integer, allocatable                :: NDend_b(:)
    Integer, allocatable                :: NDend_l(:)
    Integer, allocatable                :: Tab_ind(:)
    Integer, allocatable                :: Tab_ind1(:)
    Integer, allocatable                :: Tab_li(:)
    Integer, allocatable                :: Tags1(:)
    logical                             :: Endloop
    Integer                             :: err_io,nb,nq,i,iq,j,nb_basis,LB,inb,D
    Integer                             :: A_smol,B_smol,S,n,li_int
    Character (len=Name_len)            :: name
    Real(kind=Rk)                       :: A,B,scaleQ,Q0,d0,d2,X1,W1

    IF (debug) THEN
     Write(out_unitp,*) 'BEGINNING Read_Construct_Basis'
     Call Write_basis(Basis)
     flush(out_unitp)
    END IF


    Call Read_Basis_simple(name,nb_basis,LB,A_smol,B_smol,nb,nq,A,B,scaleQ,Q0,nio)
    Basis%Basis_name        = name
    Call string_uppercase_TO_lowercase(Basis%Basis_name)
    IF (nb_basis > 1) THEN

      Basis%NDindexq%sys_type = name
      Basis%NDindexb%sys_type = name
      Basis%NDindexl%sys_type = name

      Basis%NDindexl%NDim = nb_basis
      Basis%NDindexq%NDim = nb_basis
      Basis%NDindexb%NDim = nb_basis
      Basis%NDindexl%L        = LB


      Call string_uppercase_TO_lowercase(Basis%NDindexl%sys_type)
      Call string_uppercase_TO_lowercase(Basis%NDindexq%sys_type)
      Call string_uppercase_TO_lowercase(Basis%NDindexb%sys_type)

      SELECT CASE (Basis%Basis_name)
      CASE ('dp')

        Allocate(Basis%tab_basis(nb_basis))
        Allocate(NDend_q(nb_basis))
        Allocate(Basis%NDindexq%NDinit(nb_basis))
        Allocate(NDend_b(nb_basis))
        Allocate(Basis%NDindexb%NDinit(nb_basis))
        !iq=0
        DO i=1,nb_basis
         Call Read_Basis_simple(name,nb_basis,LB,A_smol,B_smol,nb,nq,A,B,scaleQ,Q0,nio)
         Basis%tab_basis(i)%nb_basis  = nb_basis
         Basis%tab_basis(i)%nb        = nb
         Basis%tab_basis(i)%nq        = nq
         Basis%tab_basis(i)%A         = A
         Basis%tab_basis(i)%B         = B
         Basis%tab_basis(i)%scaleQ    = scaleQ
         Basis%tab_basis(i)%Q0        = Q0
         Basis%tab_basis(i)%Basis_name= trim(adjustl(name))

         Call Construct_basis_new (Basis%tab_basis(i))
        END DO

        Basis%nb  = product(Basis%tab_basis(:)%nb)
        Basis%nq  = product(Basis%tab_basis(:)%nq)
        Basis%ndim= sum(Basis%tab_basis(:)%ndim)

        Allocate(Basis%tab_ic_to_i_dim(Basis%ndim))
        Allocate(Basis%tab_ic_to_ib(Basis%ndim))
        iq=0
        DO i=1,Size(Basis%tab_basis)
          NDend_q(i)=Basis%tab_basis(i)%nq
          NDend_b(i)=Basis%tab_basis(i)%nb
          Call Cores_ind(Basis%tab_basis(i),iq)
          Call mapp_idim(Basis,i)

        END DO

        iq=0
        DO i=1,Size(Basis%tab_basis)
          iq=iq+1
          IF(Basis%tab_basis(i)%ndim==1)THEN
            Basis%tab_ic_to_ib(iq)=i
          ELSE
            Basis%tab_ic_to_ib(iq)=i
            Basis%tab_ic_to_ib(iq+1)=i
            iq=iq+1
          END IF
        END DO

        Write(out_unitp,*)'Basis%tab_ic_to_ib(:)',Basis%tab_ic_to_ib(:)
        Write(out_unitp,*)'Basis%tab_ic_to_i_dim(:)',Basis%tab_ic_to_i_dim(:)

        Call Init_NDindex(Basis%NDindexq,NDend_q,Size(Basis%tab_basis))
        Call Init_NDindex(Basis%NDindexb,NDend_b,Size(Basis%tab_basis))

      CASE ('smolyak')

        Allocate(NDend_l(Basis%NDindexl%NDim))
        Allocate(Tab_ind(Basis%NDindexl%NDim))
        Allocate(Tab_ind1(Basis%NDindexl%NDim))
        Allocate(Tags1(Basis%NDindexl%NDim+1))

        NDend_l(:)=Basis%NDindexl%L
        Call Init_NDindex(Basis%NDindexl,NDend_l,Basis%NDindexl%NDim)
        Call Init_tab_ind(Tab_ind,Basis%NDindexl)
        I=0
        DO
          I=I+1
          CALL increase_NDindex(Tab_ind,Basis%NDindexl,Endloop)
          Basis%NDindexl%Nterm = I
          IF (Endloop) exit
        END DO
        Deallocate(Tab_ind)

        Allocate(Basis%Smolyak%DSmol(Basis%NDindexl%Nterm))

        Call Init_tab_ind(Tab_ind1,Basis%NDindexl)
        I=0
        DO
          I=I+1
          CALL increase_NDindex(Tab_ind1,Basis%NDindexl,Endloop)
          CALL Weight_D_smol (Basis%Smolyak%DSmol(I),Basis%NDindexl%L,sum(Tab_ind1),Basis%NDindexl%NDim)
          IF (Endloop) exit
        END DO
        Deallocate(Tab_ind1)

        Allocate(Basis%tab_Smolyak(Basis%NDindexl%Nterm))

        Allocate(Tab_ind(Basis%NDindexl%NDim))

        Call Init_tab_ind(Tab_ind,Basis%NDindexl)
        I=0
        DO
          I=I+1
          CALL increase_NDindex(Tab_ind,Basis%NDindexl,Endloop)

          Allocate(Basis%tab_Smolyak(I)%NDindexb%Li_smol(Basis%NDindexl%NDim))

          DO inb=1,Basis%NDindexl%NDim
            Basis%tab_Smolyak(I)%NDindexb%Li_smol(inb) = tab_ind(inb)

          END DO

          IF (Endloop) exit

        END DO
        Deallocate(Tab_ind)


        DO I=1,Basis%NDindexl%Nterm
          Allocate(Basis%tab_Smolyak(I)%tab_basis(Basis%NDindexl%NDim))
        END DO


        DO J=1, Basis%NDindexl%NDim
         Call Read_Basis_simple(name,nb_basis,LB,A_smol,B_smol,nb,nq,A,B,scaleQ,Q0,nio)

         DO I=1,Basis%NDindexl%Nterm
           Basis%tab_Smolyak(I)%Basis_name = 'dp'
           Basis%tab_Smolyak(I)%tab_basis(J)%nb_basis= nb_basis
           Basis%tab_Smolyak(I)%tab_basis(J)%A_smol  = A_smol
           Basis%tab_Smolyak(I)%tab_basis(J)%B_smol  = B_smol
           Basis%tab_Smolyak(I)%tab_basis(J)%A       = A
           Basis%tab_Smolyak(I)%tab_basis(J)%B       = B
           Basis%tab_Smolyak(I)%tab_basis(J)%scaleQ  = scaleQ
           Basis%tab_Smolyak(I)%tab_basis(J)%Q0      = Q0
           Basis%tab_Smolyak(I)%tab_basis(J)%nq = Calc_n_smol(A_smol,B_smol,Basis%tab_Smolyak(I)%NDindexb%Li_smol(J))
           Basis%tab_Smolyak(I)%tab_basis(J)%nb = Calc_n_smol(A_smol,B_smol,Basis%tab_Smolyak(I)%NDindexb%Li_smol(J))
           Basis%tab_Smolyak(I)%tab_basis(J)%Basis_name = trim(adjustl(name))

           li_int=Basis%tab_Smolyak(I)%NDindexb%Li_smol(J)
           Call Construct_basis_new (Basis%tab_Smolyak(I)%tab_basis(J),li_int)


         END DO

        END DO

        DO I=1,Basis%NDindexl%Nterm

          Basis%tab_Smolyak(I)%ndim= sum(Basis%tab_Smolyak(I)%tab_basis(:)%ndim)
          Basis%tab_Smolyak(I)%nb = product(Basis%tab_Smolyak(I)%tab_basis(:)%nb)
          Basis%tab_Smolyak(I)%nq = product(Basis%tab_Smolyak(I)%tab_basis(:)%nq)
          Basis%tab_Smolyak(I)%NDindexb%sys_type = 'dp'
          Basis%tab_Smolyak(I)%NDindexq%sys_type = 'dp'

          Allocate(NDend_q(Basis%NDindexl%NDim))
          Allocate(NDend_b(Basis%NDindexl%NDim))

          DO J=1,Basis%NDindexl%NDim
            NDend_q(J) = Basis%tab_Smolyak(I)%tab_basis(J)%nq
            NDend_b(J) = Basis%tab_Smolyak(I)%tab_basis(J)%nb

          END DO


          Allocate(Basis%tab_Smolyak(I)%tab_ic_to_i_dim(Basis%tab_Smolyak(I)%ndim))
          Allocate(Basis%tab_Smolyak(I)%tab_ic_to_ib(Basis%tab_Smolyak(I)%ndim))

          Iq=0
          DO n=1,Size(Basis%tab_Smolyak(I)%tab_basis)
              Call Cores_ind(Basis%tab_Smolyak(I)%tab_basis(n),iq)
              Call mapp_idim(Basis%tab_Smolyak(I),n)
          END DO

          Iq = 0
          DO n=1,Size(Basis%tab_Smolyak(I)%tab_basis)
            Iq = Iq+1
            IF(Basis%tab_Smolyak(I)%tab_basis(n)%ndim==1)THEN
              Basis%tab_Smolyak(I)%tab_ic_to_ib(iq)=n
            ELSE
              Basis%tab_Smolyak(I)%tab_ic_to_ib(iq)=n
              Basis%tab_Smolyak(I)%tab_ic_to_ib(iq+1)=n
              Iq=Iq+1
            END IF
          END DO

          Call Init_NDindex(Basis%tab_Smolyak(I)%NDindexq,NDend_q,Basis%NDindexl%NDim)
          Call Init_NDindex(Basis%tab_Smolyak(I)%NDindexb,NDend_b,Basis%NDindexl%NDim)

          Deallocate(NDend_q)
          Deallocate(NDend_b)

        END DO

        Call  Mapping_2G (Basis )
        CALL test_smolyak(Basis)
      CASE default
        STOP 'ERROR in Read_Basis: no default basis.'
      END SELECT

    ELSE
      Basis%nb_basis  = nb_basis
      Basis%nb        = nb
      Basis%nq        = nq
      Basis%A         = A
      Basis%B         = B
      Basis%scaleQ    = scaleQ
      Basis%Q0        = Q0
      Call Construct_basis_new (Basis)
    END IF
    IF (debug) THEN
     Write(out_unitp,*) 'END Read_Construct_Basis'
     Call Write_basis(Basis)
     flush(out_unitp)
    END IF
  END SUBROUTINE read_construct_basis


  SUBROUTINE Calc_nbi_li(Basis)
    USE UtilLib_m
    IMPLICIT NONE
     !Logical,   parameter          :: debug = .true.
     Logical,    parameter          :: debug = .false.
     TYPE(Basis_t),intent(inout)    :: Basis
     integer                        :: i,j,N

     IF (debug) THEN
      Write(out_unitp,*) 'Calc_nbi_li'
      flush(out_unitp)
     END IF
     !Termsmo|l1| l2 |n1|n2|nb|
    !    1   |0 | 0  |1 |1 |1 |
    !    2   |1 | 0  |9 |1 |9 |
    !    3   |2 | 0  |25|1 |25|
    !    4   |0 | 1  |1 |3 |3 |
    !    5   |1 | 1  |9 |3 |27|
    !    6   |0 | 2  |1 |5 |5 |
     Allocate(Basis%Nbi_li(0:Basis%NDindexl%L,Basis%NDindexl%NDim))

      DO n = 1,Basis%NDindexl%Nterm
        DO j = 1,Basis%NDindexl%NDim!  j=1 et 2
          i = Basis%tab_Smolyak(n)%NDindexb%Li_smol(j)
          Basis%Nbi_li(i,j) = Basis%tab_Smolyak(n)%tab_basis(J)%nb

        END DO
      END DO

     IF (debug) THEN
       Write(out_unitp,*)'Basis%Nbi_li='
       DO i=0,Basis%NDindexl%L
        Write(out_unitp,*) Basis%Nbi_li(i,:)
       END DO
      Write(out_unitp,*) 'Calc_nbi_li'
      flush(out_unitp)
     END IF
  END SUBROUTINE Calc_nbi_li


  SUBROUTINE Construct_basis_new (Basis,l)
   USE UtilLib_m
   IMPLICIT NONE
     !Logical,            parameter          :: debug = .true.
     Logical,             parameter          :: debug = .false.
     TYPE(Basis_t),  intent(inout)           :: Basis
     Integer,       intent(in), optional     :: l
     Integer                                 :: li


     IF (debug) THEN
      Write(out_unitp,*) 'Begining Construct_basis_new'
      Call write_basis(Basis)
      flush(out_unitp)
     END IF

     IF (Present(l)) THEN
       li=l
     ELSE
       li=-2
     END IF

     Call string_uppercase_TO_lowercase(Basis%Basis_name)
     SELECT CASE (Basis%Basis_name)
      CASE ('boxab')
       Basis%ndim=1
       Call Construct_Basis_Sin(Basis)
       Basis%Q0      = Basis%A
       Basis%scaleQ  = pi/(Basis%B-Basis%A)
       Basis%ndim=1

      CASE ('herm','ho')
       Basis%ndim=1
       Call Construct_Basis_Ho(Basis)

      CASE ('dim_2')
       Basis%ndim=2
       Call Construct_Basis_2B(Basis,li)

      CASE ('four','fourier')
       Call Construct_Basis_Fourier(Basis)
      CASE default
       STOP 'ERROR in Read_Basis: no default basis.'
     END SELECT

     Call Scale_Basis(Basis)
     Call Calc_dngg_grid(Basis)
     Call CheckOrtho_Basis(Basis,nderiv=2)

     IF (debug) THEN
      Call write_basis(Basis)
      Write(out_unitp,*) 'END Construct_basis_new'
      flush(out_unitp)
     END IF

  END SUBROUTINE Construct_basis_new

  SUBROUTINE Construct_Basis_2B(Basis,l) ! 2basis
   USE UtilLib_m
   IMPLICIT NONE
     TYPE(Basis_t),  intent(inout)       :: Basis
     TYPE(Basis_t),allocatable           :: Basis_2(:)
     Real(kind=Rk)                       :: dx
     Integer                             :: ib,iq,nb,nq
     Integer , intent(in)                :: l
     Integer                             :: ib1,ib2,iq1,iq2,i
     !Logical,         parameter         :: debug = .true.
     Logical, parameter                  :: debug = .false.

     IF (debug) THEN

      Write(out_unitp,*) 'Begining Construct_Basis_2B'
      Call write_basis(Basis)
      flush(out_unitp)
    END IF

    IF (l/=-2) THEN
       Basis%nq = (Calc_n_smol(Basis%A_smol,Basis%B_smol,l))**2
       Basis%nb = (Calc_n_smol(Basis%A_smol,Basis%B_smol,l))**2
       nq=Basis%nq
       nb=Basis%nb
    ELSE
        nq=Basis%nq
        nb=Basis%nb
    END IF

    Allocate(Basis_2(2))

    Do i=1,2
       Basis_2(i)%nq        = Int(sqrt(Real(Basis%nq)))

       Basis_2(i)%nb        = int(sqrt(Real(Basis%nb)))
       Basis_2(i)%nb_basis  = Basis%nb_basis
       Basis_2(i)%Ndim      = 1
       Basis_2(i)%A         = Basis%A
       Basis_2(i)%B         = Basis%B
       Basis_2(i)%scaleQ    = Basis%scaleQ
       Basis_2(i)%Q0        = Basis%Q0
       Basis_2(i)%Basis_name= 'herm'
       Call Construct_basis_new (Basis_2(i))

     END DO


     Allocate(Basis%W(nq))
     Allocate(Basis%x(nq,2))
     Allocate(Basis%d0gb(nq,nb))
     Allocate(Basis%d1gb(nq,nb,2))
     Allocate(Basis%d2gb(nq,nb,2,2))

     ib1=1
     ib2=0
     DO Ib=1,nb
       IF (ib2 == Basis_2(2)%nb) THEN
         ib1 = ib1 + 1
         ib2 = 1
       ELSE
         ib2 = ib2 + 1
       END IF
         iq1=1
         iq2=0
         DO Iq=1,nq
          IF (iq2 == Basis_2(2)%nq) THEN
           iq1 = iq1 + 1
           iq2 = 1
          ELSE
           iq2 = iq2 + 1
          END IF


          Basis%W(Iq) = Basis_2(1)%w(iq1)*Basis_2(2)%w(iq2)

          Basis%x(iq,1) = Basis_2(1)%x(iq1,1)
          Basis%x(iq,2) = Basis_2(2)%x(iq2,1)

          Basis%d0gb(iq,ib) = Basis_2(1)%d0gb(iq1,ib1)*Basis_2(2)%d0gb(iq2,ib2)
          Basis%d1gb(iq,ib,1) = Basis_2(1)%d1gb(iq1,ib1,1)*Basis_2(2)%d0gb(iq2,ib2)
          Basis%d1gb(iq,ib,2) = Basis_2(1)%d0gb(iq1,ib1)*Basis_2(2)%d1gb(iq2,ib2,1)

          Basis%d2gb(iq,ib,1,1) = Basis_2(1)%d2gb(iq1,ib1,1,1)*Basis_2(2)%d0gb(iq2,ib2)
          Basis%d2gb(iq,ib,1,2) = Basis_2(1)%d1gb(iq1,ib1,1)*Basis_2(2)%d1gb(iq2,ib2,1)
          Basis%d2gb(iq,ib,2,2) = Basis_2(1)%d0gb(iq1,ib1)*Basis_2(2)%d2gb(iq2,ib2,1,1)

         END DO
       END DO

       IF (debug) THEN
        Call write_basis(Basis)
        Write(out_unitp,*) 'END Construct_Basis_2B'
        flush(out_unitp)
      END IF
  END SUBROUTINE Construct_Basis_2B

  SUBROUTINE Check_index(I_smol,n,Tab_ind,Basis)
   IMPLICIT NONE
   TYPE(Basis_t), intent(inout)  :: Basis
   Integer,      intent(in)      :: Tab_ind(:),n,I_smol
   Integer                       :: I_compt

   DO I_compt = 1,Basis%NB
      IF( ALL(Tab_ind(1:Basis%NDindexl%NDim) == Basis%tab_in_Smolyak(I_compt)%Tags1(1:Basis%NDindexl%NDim)))THEN
        Basis%tab_Smolyak(I_smol)%Ind_map(n) = I_compt

        exit
      END IF
   END DO
  END SUBROUTINE Check_index


  SUBROUTINE Mapping_2G (Basis )
    IMPLICIT NONE
    !Logical,     parameter      :: debug = .true.
    Logical,     parameter       :: debug = .false.
    TYPE(Basis_t),intent(inout)  :: Basis
    Integer                      :: I_smol,n,i,I_smol1,compt
    Integer, allocatable         :: Tab_ind(:)
    Integer, allocatable         :: Tab_li(:)
    Integer, allocatable         :: Ndend(:,:)
    Integer, allocatable         :: Ndbeging(:)
    Integer                      :: Nbinter,I_smol2,I_compt
    Integer, allocatable         :: Nbinter1(:)
    Integer, allocatable         :: Tags(:)
    Integer, allocatable         :: Tags1(:)
    Logical                      :: Endloop

    Call Calc_nbi_li(Basis)

    Allocate(Tab_li(Basis%NDindexl%NDim))
    Allocate(Ndbeging(Basis%NDindexl%NDim))
    Allocate(Nbinter1(Basis%NDindexl%NDim))
    Allocate(Ndend(Basis%NDindexl%Nterm,Basis%NDindexl%NDim))

    Call Init_tab_ind(Tab_li,Basis%NDindexl)
    I_smol = 0
    Compt  = 0

    DO
      I_smol = I_smol+1

      CALL increase_NDindex(Tab_li,Basis%NDindexl,Endloop)

      Allocate(Tags1(Basis%NDindexl%NDim))

      DO i=1,Basis%NDindexl%NDim
        Ndend(I_smol,i) = Basis%tab_Smolyak(I_smol)%tab_basis(i)%nb
      END DO


      DO i = 1,Basis%NDindexl%NDim
        IF(Tab_li(i)==0)THEN
          Ndbeging(i) =   Basis%Nbi_li(Tab_li(i),i)
          Nbinter1(i) =   Basis%Nbi_li(Tab_li(i),i)
        ELSE
          Ndbeging(i) = Basis%Nbi_li(Tab_li(i)-1,i)+1
          Nbinter1(i) = Basis%Nbi_li(Tab_li(i),i)-Basis%Nbi_li(Tab_li(i)-1,i)
        END IF
      END DO

      Nbinter  = Product(Nbinter1(:))

      Tags1(1:Basis%NDindexl%NDim) = Ndbeging(:)
      Tags1(1) = Ndbeging(1)-1

      DO  n = 1,Nbinter
        Tags1(1) = Tags1(1)+1
        DO i = 1,Basis%NDindexl%NDim-1

          IF(Tags1(i) > Ndend(I_smol,i)) THEN
           Tags1(i+1) = Tags1(i+1)+1
           Tags1(i)   = Ndbeging(i)
          END IF

        END DO

        Compt = Compt+1

      END DO
      Deallocate(Tags1)
      IF (Endloop) exit
    END DO


    Allocate(Basis%tab_in_Smolyak(Compt))

    Call Init_tab_ind(Tab_li,Basis%NDindexl)

    I_smol2 = 0
    Compt   = 0

    DO
      I_smol2 = I_smol2 + 1

      CALL increase_NDindex(Tab_li,Basis%NDindexl,Endloop)

      Allocate(Tags(Basis%NDindexl%NDim+1))


      DO i = 1,Basis%NDindexl%NDim
        IF(Tab_li(i)==0)THEN
         Ndbeging(i) =   Basis%Nbi_li(Tab_li(i),i)
         Nbinter1(i) =   Basis%Nbi_li(Tab_li(i),i)
        ELSE
         Ndbeging(i) = Basis%Nbi_li(Tab_li(i)-1,i)+1
         Nbinter1(i) = Basis%Nbi_li(Tab_li(i),i)-Basis%Nbi_li(Tab_li(i)-1,i)
        END IF
      END DO

      Nbinter  = Product(Nbinter1(:))
      Tags(1:Basis%NDindexl%NDim) = Ndbeging(:)
      Tags(1) = Ndbeging(1)-1

      DO  n = 1,Nbinter
        Tags(1) = Tags(1)+1
        DO i = 1,Basis%NDindexl%NDim-1
          IF(Tags(i)  >  Ndend(I_smol2,i)) THEN
           Tags(i+1) =  Tags(i+1)+1
           Tags(i)   =  Ndbeging(i)
          END IF
        END DO
        Compt = Compt+1
        Allocate(Basis%tab_in_Smolyak(Compt)%Tags1(Basis%NDindexl%NDim+1))
        Tags(Basis%NDindexl%NDim+1)        = Compt
        Basis%tab_in_Smolyak(Compt)%Tags1  = Tags
      END DO
      Deallocate(Tags)
      IF (Endloop) exit
    END DO

    Basis%NB = compt
    Write(*,*) 'Nb_compate=',Basis%NB
    Allocate(Tab_ind(Basis%NDindexl%NDim))

    DO I_smol=1,Basis%NDindexl%Nterm
      Allocate(Basis%tab_Smolyak(I_smol)%Ind_map(Basis%tab_Smolyak(I_smol)%nb))
      Tab_ind(:) = 1
      Tab_ind(1) = 0
      DO  n = 1,Basis%tab_Smolyak(I_smol)%nb
        CALL increase_NDindex_simple(Tab_ind,Ndend(I_smol,:),Basis%NDindexl)
        CALL Check_index(I_smol,n,Tab_ind,Basis)
      END DO
    END DO

  END SUBROUTINE Mapping_2G


  SUBROUTINE test_smolyak(Basis)
   USE UtilLib_m
   IMPLICIT NONE
    TYPE(Basis_t), intent(inout)  :: Basis
    Integer                       :: I_smol,n,i,I_smol1,I_compt
    Real(kind=Rk), allocatable    :: Vsort(:),V_entre(:),Vb(:),Vg(:)
    !Logical,parameter             :: debug = .true.
    Logical, parameter            :: debug = .false.

    IF (debug) THEN
     Write(out_unitp,*)'BEGINNING test_smolyak'
     Write(out_unitp,*) 'Basis%NB',Basis%NB
     flush(out_unitp)
    END IF

    Allocate(V_entre(Basis%NB ))
    Allocate(Vsort(Basis%NB ))

    V_entre(:) = Three
    V_entre(1) = One
    Vsort(:)=ZERO

    DO I_smol = 1,Basis%NDindexl%Nterm

      Allocate(Vb(Basis%tab_Smolyak(I_smol)%nb))
      Allocate(Vg(Basis%tab_Smolyak(I_smol)%nq))

      DO n = 1,Basis%tab_Smolyak(I_smol)%nb
        I_compt = Basis%Tab_Smolyak(I_smol)%Ind_map(n)
        Vb(n) = V_entre(I_compt)
      END DO

      Call BasisTOGrid_Basis_rapide1(Vg,Vb, Basis%tab_Smolyak(I_smol))
      Call GridTOBasis_Basis_rapide1(Vb,Vg, Basis%tab_Smolyak(I_smol))

      DO n = 1,Basis%tab_Smolyak(I_smol)%nb
        I_compt = Basis%Tab_Smolyak(I_smol)%Ind_map(n)
        Vsort(I_compt)= Vsort(I_compt)+ Basis%Smolyak%DSmol(I_smol)*Vb(n)
      END DO
      Deallocate(Vb)
      Deallocate(Vg)

    END DO

    IF (debug) THEN
     Write(out_unitp,*)'end test_smolyak'

     flush(out_unitp)
    END IF
  END SUBROUTINE test_smolyak

  SUBROUTINE Cores_ind(Basis,iq)
    IMPLICIT NONE
    TYPE(Basis_t),intent(inout)  :: Basis
    Integer,intent(inout)        :: iq
    Integer                      :: icompt
    !Logical,parameter            :: debug = .true.
    Logical, parameter        :: debug = .false.

     IF (debug) THEN
      Write(out_unitp,*)'BEGINNING Cores_ind'
      Write(out_unitp,*)"ndim=",Basis%ndim,"iq=",iq
      flush(out_unitp)
     END IF

     Allocate(Basis%ind_iq(Basis%ndim))
     Do icompt=1,Basis%ndim
        Basis%ind_iq(icompt)= iq+icompt
     END DO

     Iq = Basis%ind_iq(Basis%ndim)

     IF (debug) THEN
       Write(out_unitp,*) 'END Cores_ind'
       Write(out_unitp,*) "iq=",iq
        DO icompt=1,Basis%ndim
         Write(out_unitp,*) "ind_iq=",Basis%ind_iq(icompt)
        END   DO
       flush(out_unitp)
     END IF
  END SUBROUTINE Cores_ind


  SUBROUTINE mapp_idim(Basis,i)
    IMPLICIT NONE
    TYPE(Basis_t),intent(inout)  :: Basis
    Integer,intent(in)           :: i
    Integer                      :: icompt
    !Logical,parameter            :: debug = .true.
    Logical, parameter        :: debug = .false.

     IF (debug) THEN
      Write(out_unitp,*)'BEGINNING mapp_idim'
      Write(out_unitp,*)"ndim=",Basis%tab_basis(i)%ndim
      flush(out_unitp)
     END IF

     Do icompt=1,Basis%tab_basis(i)%ndim
        Basis%tab_ic_to_i_dim(Basis%tab_basis(i)%ind_iq(icompt))=icompt
     END DO
     !Write(out_unitp,*) 'i=',i
    IF (debug) THEN
    Write(out_unitp,*) '----------------------------------------------------------'
     Write(out_unitp,*) 'END mapp_idim'
     DO icompt=1,Basis%tab_basis(i)%ndim
       !Write(out_unitp,*) "ind_iq=",Basis%ind_iq(icompt)
        Write(out_unitp,*) "tab_ic_to_i_dim=",&
        Basis%tab_ic_to_i_dim(Basis%tab_basis(i)%ind_iq(icompt))
     END   DO
     flush(out_unitp)
    END IF
  END SUBROUTINE mapp_idim



 SUBROUTINE Weight_D_smol ( D,L,Som_l,NDIM)
  IMPLICIT NONE
  Integer,   intent(in)       :: NDIM
  Integer,   intent(in)       :: Som_l
  Integer,   intent(in)       :: L
  Integer,   intent(out)      :: D

  D = ((-1)**(L-Som_l))*Binomial(L-Som_l,NDIM-1)

 END SUBROUTINE Weight_D_smol


 SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
 USE UtilLib_m
 IMPLICIT NONE
   TYPE(Basis_t),       intent(inout)  :: Basis
   Real(kind=Rk)                       :: dx
   Integer                             :: ib,iq,nb,nq

   nb = Basis%nb
   nq = Basis%nq
   dx = pi/nq
   Allocate(Basis%W(nq))
   Allocate(Basis%x(nq,1))
   Allocate(Basis%d0gb(nq,nb))
   Allocate(Basis%d1gb(nq,nb,1))
   Allocate(Basis%d2gb(nq,nb,1,1))

   Basis%x(:,1) = [(dx*(iq-HALF),iq=1,nq)]
   Basis%w(:) = [(dx,iq=1,nq)]


   DO ib=1,nb
     Basis%d0gb(:,ib)     =          sin(Basis%x(:,1)*ib) / sqrt(pi*HALF)
     Basis%d1gb(:,ib,1)   =  ib    * cos(Basis%x(:,1)*ib) / sqrt(pi*HALF)
     Basis%d2gb(:,ib,1,1) = -ib**2 * Basis%d0gb(:,ib)
   END DO

   IF (nb == nq) THEN
     Basis%d0gb(:,nb)      = Basis%d0gb(:,nb)      / sqrt(TWO)
     Basis%d1gb(:,nb,:)    = Basis%d1gb(:,nb,:)    / sqrt(TWO)
     Basis%d2gb(:,nb,:,:)  = Basis%d2gb(:,nb,:,:)  / sqrt(TWO)
   END IF
 END SUBROUTINE Construct_Basis_Sin

 SUBROUTINE Construct_Basis_Fourier(Basis)
  USE UtilLib_m
  IMPLICIT NONE
   TYPE(Basis_t),  intent(inout)  :: Basis
   Real(kind=Rk)                  :: dx
   Integer                        :: ib,iq,nb,nq,K

   nb = Basis%nb
   nq = Basis%nq
   dx = TWO*PI/nq
   Allocate(Basis%W(nq))
   Allocate(Basis%x(nq,1))
   Allocate(Basis%d0gb(nq,nb))
   Allocate(Basis%d1gb(nq,nb,1))
   Allocate(Basis%d2gb(nq,nb,1,1))

   Basis%x(:,1) = [(dx*iq-dx/2 , iq = 1,nq)]
   Basis%w(:)   = [(dx,iq = 1,nq)]

   DO ib = 1, nb
    K= int(ib/2)
    IF (ib == 1) THEN
      Basis%d0gb(:, ib)        =  ONE/sqrt(TWO*PI)
      Basis%d1gb(:, ib, 1)     =  ZERO
      Basis%d2gb(:, ib, 1, 1)  =  ZERO
    ELSE IF (modulo(ib, 2) == 0) THEN
      Basis%d0gb(:, ib)        =  sin(Basis%x(:,1)*K)/sqrt(PI)
      Basis%d1gb(:, ib, 1)     =  K*cos(Basis%x(:,1)*K)/sqrt(PI)
      Basis%d2gb(:, ib, 1, 1)  = -K**2*Basis%d0gb(:, ib)
    ELSE
      Basis%d0gb(:, ib)        =  cos(Basis%x(:,1)*K)/sqrt(PI)
      Basis%d1gb(:, ib, 1)     = -K*sin(Basis%x(:,1)*K)/sqrt(PI)
      Basis%d2gb(:, ib, 1, 1)  = -K**2*Basis%d0gb(:, ib)
    END IF
   END DO

   IF (nb == nq .and. modulo(Basis%nb,2) == 0) THEN
     Basis%d0gb(:,nb)          =  Basis%d0gb(:,nb)      / sqrt(TWO)
     Basis%d1gb(:,nb,:)        =  Basis%d1gb(:,nb,:)    / sqrt(TWO)
     Basis%d2gb(:,nb,:,:)      =  Basis%d2gb(:,nb,:,:)  / sqrt(TWO)
   END IF

 END SUBROUTINE Construct_Basis_Fourier


 SUBROUTINE Construct_Basis_Ho(Basis) ! HO :
 USE UtilLib_m
  TYPE(Basis_t), intent(inout)  :: Basis
  Integer                       :: iq,ib

!  Allocate(Basis%x(Basis%nq))
  Allocate(Basis%x(Basis%nq,1))
  Allocate(Basis%w(Basis%nq))

  Call hercom(Basis%nq, Basis%x(:,1), Basis%w(:))

  Allocate(Basis%d0gb(Basis%nq,Basis%nb))
  Allocate(Basis%d1gb(Basis%nq,Basis%nb,1))
  Allocate(Basis%d2gb(Basis%nq,Basis%nb,1,1))

  DO iq = 1, Basis%nq
    DO ib = 1, Basis%nb
      Call Construct_Basis_poly_Hermite_exp(Basis%x(iq,1),Basis%d0gb(iq,ib),&
      Basis%d1gb(iq,ib,1),Basis%d2gb(iq,ib,1,1), ib-1,.TRUE.)
    END DO
  END DO
 END SUBROUTINE Construct_Basis_Ho


 SUBROUTINE CheckOrtho_Basis(Basis,nderiv)
 USE UtilLib_m
   !Logical,                 parameter   :: debug = .true.
   Logical,                parameter   ::debug = .false.
   TYPE(Basis_t),           intent(in)   :: Basis
   Integer,                 intent(in)   :: nderiv
   Integer                               :: ib
   Real(kind=Rk), ALLOCATABLE            :: S(:,:)
   Real(kind=Rk), ALLOCATABLE            :: d0bgw(:,:)
   Real(kind=Rk)                         :: Sii,Sij

   IF (Basis_IS_Allocated(Basis)) THEN
    d0bgw = transpose(Basis%d0gb)
    DO ib=1,Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
    END DO
    S = matmul(d0bgw,Basis%d0gb)
    !IF (nderiv > -1) Call Write_RMat(S,out_unitp,5,name_info='S')
    Sii = ZERO
    Sij = ZERO
    DO ib=1,Basis%nb
      IF (abs(S(ib,ib)-ONE) > Sii) Sii = abs(S(ib,ib)-ONE)
      S(ib,ib) = ZERO
    END DO

    Sij = maxval(S)
  !  write(out_unitp,*) 'Sii,Sij',Sii,Sij

    IF (nderiv > 0) THEN
      write(out_unitp,*)
      S = matmul(d0bgw,Basis%d1gb(:,:,1))
      !Call Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>',Rformat='e13.4')
    !  Call Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
    END IF

    IF (nderiv > 1) THEN
      write(out_unitp,*)
      S = matmul(d0bgw,Basis%d2gb(:,:,1,1))
      !Call Write_RMat(S,out_unitp,5,name_info='<d0b|d2b>',Rformat='e13.4')
    !  Call Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
    END IF
  ELSE
    write(out_unitp,*) ' WARNNING in CheckOrtho_Basis'
    write(out_unitp,*) ' the basis is not Allocated.'
  END IF
  END SUBROUTINE CheckOrtho_Basis

  SUBROUTINE BasisTOGrid_Basis(G,B,Basis)
  USE UtilLib_m
  USE NDindex_m
    !Logical,           parameter    :: debug = .true.
    Logical,           parameter     :: debug = .false.
    TYPE(Basis_t),     intent(in)    :: Basis
    Real(kind=Rk),     intent(in)    :: B(:) !Vector on base,
    Real(kind=Rk),     intent(inout) :: G(:) !Vector on the grid, out
    Real(kind=Rk)                    :: W
    Logical                          :: Endloop_q
    Logical                          :: Endloop_b
    Integer,         allocatable     :: tab_iq(:)
    Integer,         allocatable     :: tab_ib(:)
    Integer                          :: ib,iq,nq,nb,inb
    Integer                          :: jb,jb1,jb2

    IF(debug)THEN
     write(out_unitp,*) 'BEGINNING BasisTOGrid_Basis'
     write(out_unitp,*) 'intent(in) :: B(:)',B
     Call Write_Basis(Basis)
     flush(out_unitp)
    END IF

    IF(.NOT. Basis_IS_Allocated(Basis))THEN
     write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
     write(out_unitp,*) " the basis is not Allocated."
     STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
    END IF

    IF(size(B) /= Basis%nb)THEN
     write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
     write(out_unitp,*) ' the size of B is different from nb.'
     write(out_unitp,*) ' size(B), Basis%nb',size(B),Basis%nb
     STOP 'ERROR in BasisTOGrid_Basis: wrong B size.'
    END IF

    IF(size(G) /= Basis%nq)THEN
     write(out_unitp,*) ' ERROR in GridTOBasis_Basis'
     write(out_unitp,*) ' the size of G is different from nq.'
     write(out_unitp,*) ' size(G), Basis%nq',size(G),Basis%nq
     STOP 'ERROR in BasisTOGrid_Basis: wrong G size..'
    END IF

    IF(Allocated(Basis%tab_basis)) THEN
     Allocate(Tab_ib(size(Basis%tab_basis)))
     Allocate(Tab_iq(size(Basis%tab_basis)))
     Call Init_tab_ind(Tab_iq,Basis%NDindexq)
     Iq=0
     DO
      Iq=Iq+1
      Call increase_NDindex(Tab_iq,Basis%NDindexq,Endloop_q)
      IF (Endloop_q) exit
      G(iq)=ZERO
      Call Init_tab_ind(Tab_ib,Basis%NDindexb)
      Ib=0
      DO
       Ib=Ib+1
       Call increase_NDindex(Tab_ib,Basis%NDindexb,Endloop_b)
       IF (Endloop_b) exit
       W = ONE
       DO inb=1,size(Basis%tab_basis)
        W = W * Basis%tab_basis(inb)%d0gb(tab_iq(inb),tab_ib(inb))
       END DO
       G(iq) = G(iq) + W * B(ib)
      END DO
     END DO
     Deallocate(Tab_ib)
     Deallocate(Tab_iq)
    ELSE
     DO iq = 1,Basis%nq
      G(iq) = ZERO
      DO ib = 1,Basis%nb
       G(iq) = G(iq) + Basis%d0gb(iq,ib) * B(ib)
      END DO
     END DO
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'intent(OUTIN) :: G(:)',G
      write(out_unitp,*) 'END BasisTOGrid_Basis'
      flush(out_unitp)
    END IF
  END SUBROUTINE BasisTOGrid_Basis



  SUBROUTINE  GridTOBasis_1D(BB,GG,Basis)
   USE UtilLib_m
   TYPE(Basis_t)    , intent(in),target     :: Basis
   Real (kind=Rk), intent(inout)            :: BB(:,:,:)
   Real (kind=Rk), intent(in)               :: GG(:,:,:)
   real(kind=Rk), ALLOCATABLE               :: d0bgw(:,:)
   Logical          , parameter             :: debug = .true.
   Integer                                  :: i1,i3,ib

    IF (debug) THEN
      flush(out_unitp)
    END IF

    d0bgw = transpose(Basis%d0gb)
    DO ib=1,Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
    END DO

    BB(:,:,:) = ZERO
    DO i3=1,ubound(GG,dim=3)
    DO i1=1,ubound(GG,dim=1)

      BB(i1,:,i3) =  matmul( d0bgw  ,GG(i1,:,i3))

    END DO
    END DO

    IF (debug) THEN
      flush(out_unitp)
    END IF
  END SUBROUTINE  GridTOBasis_1D


  SUBROUTINE BasisTOGrid_1D(GB,BB,Basis)
   USE UtilLib_m
    TYPE(Basis_t)    , intent(in),target     :: Basis
    Real (kind=Rk), intent(inout)            :: GB(:,:,:)
    Real (kind=Rk), intent(in)               :: BB(:,:,:)
    Logical          , parameter             :: debug = .true.
    Integer                                  :: i1,i3,iq,ib

      IF (debug) THEN
        flush(out_unitp)
      END IF

      DO i3 = 1,ubound(BB, dim=3)
      DO i1 = 1,ubound(BB, dim=1)

        GB(i1,:,i3) =  matmul( Basis%d0gb , BB(i1,:,i3))

      END DO
      END DO


      IF (debug) THEN
      	flush(out_unitp)
      END IF
  END SUBROUTINE  BasisTOGrid_1D


  SUBROUTINE Calc_indice0( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
    TYPE(Basis_t),     intent(in),target  :: Basis
    integer,intent(inout) , allocatable   :: Ib1(:),Ib2(:),Iq3(:),Iq1(:),Iq2(:),Ib3(:)
    integer,intent(in)                    :: Ndim
    integer                               :: inb

    allocate(Ib3(Ndim))
    allocate(Ib2(Ndim))
    allocate(Ib1(Ndim))

    allocate(Iq3(Ndim))
    allocate(Iq2(Ndim))
    allocate(Iq1(Ndim))


     DO inb = 1, Ndim
       IF (inb == 1)THEN
        Iq1(1)  =  1
        Ib1(1)  =  1
        Iq2(1)  =  Basis%tab_basis(1)%nq
        Ib2(1)  =  Basis%tab_basis(1)%nb
        Iq3(1)  =  Product(Basis%tab_basis(2:Ndim)%nq)
        Ib3(1)  =  Product(Basis%tab_basis(2:Ndim)%nb)
       ELSE IF (inb == Ndim) THEN
        Iq1(inb) =  Product(Basis%tab_basis(1:Ndim-1)%nq)
        Ib1(inb) =  Product(Basis%tab_basis(1:Ndim-1)%nb)
        Ib2(inb) =  Basis%tab_basis(Ndim)%nb
        Iq2(inb) =  Basis%tab_basis(Ndim)%nq
        Ib3(inb) =  1
        Iq3(inb) =  1
       ELSE
        Iq1(inb) =  Product(Basis%tab_basis(1:inb-1)%nq)
        Ib1(inb) =  Product(Basis%tab_basis(1:inb-1)%nb)
        Iq2(inb) =  Basis%tab_basis(inb)%nq
        Ib2(inb) =  Basis%tab_basis(inb)%nb
        Ib3(inb) =  Product(Basis%tab_basis(inb+1:Ndim)%nb)
        Iq3(inb) =  Product(Basis%tab_basis(inb+1:Ndim)%nq)
       END IF
     END DO

  END SUBROUTINE Calc_indice0

  SUBROUTINE BasisTOGrid_Basis_rapide(G,B,Basis)
  USE UtilLib_m
  USE NDindex_m
    !Logical,           parameter          :: debug = .true.
    Logical,          parameter             :: debug = .false.
    TYPE(Basis_t),  intent(in)              :: Basis
    Real (kind=Rk),  intent(in),target      :: B(:) !Vector on base,
    Real (kind=Rk),  intent(inout),target   :: G(:) !Vector on the grid, out
    Real (kind=Rk), pointer                 :: BBG(:,:,:)
    Real (kind=Rk), pointer                 :: BBB(:,:,:)
    Real (kind=Rk),allocatable,target       :: GB1(:),GB2(:)
    integer , allocatable                   :: Ib3(:),Iq1(:),Iq2(:),ib1(:),ib2(:),iq3(:)
    Real (kind=Rk), pointer                 :: GBB(:,:,:)
    Real (kind=Rk) ,allocatable,target      :: GBB1(:),GGB2(:)

    Integer                                 :: ib,iq,nq,nb,inb,Ndim
    Integer                                 :: jb,jb1,jb2

    IF(debug)THEN
     write(out_unitp,*) 'BEGINNING BasisTOGrid_Basis'
     write(out_unitp,*) 'intent(in) :: B(:)',B
     Call Write_Basis(Basis)
     flush(out_unitp)
    END IF

    IF(.NOT. Basis_IS_Allocated(Basis))THEN
     write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
     write(out_unitp,*) " the basis is not Allocated."
     STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
    END IF

    IF(size(B) /= Basis%nb)THEN
     write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
     write(out_unitp,*) ' the size of B is different from nb.'
     write(out_unitp,*) ' size(B), Basis%nb',size(B),Basis%nb
     STOP 'ERROR in BasisTOGrid_Basis: wrong B size.'
    END IF

    IF(size(G) /= Basis%nq)THEN
     write(out_unitp,*) ' ERROR in GridTOBasis_Basis'
     write(out_unitp,*) ' the size of G is different from nq.'
     write(out_unitp,*) ' size(G), Basis%nq',size(G),Basis%nq
     STOP 'ERROR in BasisTOGrid_Basis: wrong G size..'
    END IF

    IF(Allocated(Basis%tab_basis)) THEN

      Allocate(GBB1(Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb*Basis%tab_basis(3)%nb))


      BBB(1:1,1:Basis%tab_basis(1)%nb,1:Basis%tab_basis(2)%nb*Basis%tab_basis(3)%nb) => B
      GBB(1:1,1:Basis%tab_basis(1)%nq,1:Basis%tab_basis(2)%nb*Basis%tab_basis(3)%nb) => GBB1

      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(1))
!

      Allocate(GGB2(Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nq*Basis%tab_basis(3)%nb))

      BBB(1:Basis%tab_basis(1)%nq,1:Basis%tab_basis(2)%nb,1:Basis%tab_basis(3)%nb) =>  GBB1
      GBB(1:Basis%tab_basis(1)%nq,1:Basis%tab_basis(2)%nq,1:Basis%tab_basis(3)%nb) =>  GGB2

      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(2))

      Deallocate(GBB1)


      BBB(1:Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nq,1:Basis%tab_basis(3)%nb,1:1)  => GGB2
      GBB(1:Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nq,1:Basis%tab_basis(3)%nq,1:1 ) => G

      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(3))


    ELSE
      DO iq = 1,Basis%nq
       G(iq) = ZERO
       DO ib = 1,Basis%nb
        G(iq) = G(iq) + Basis%d0gb(iq,ib) * B(ib)
       END DO
      END DO
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'intent(OUTIN) :: G(:)',G
      write(out_unitp,*) 'END BasisTOGrid_Basis_rapide'
      flush(out_unitp)
    END IF
  END SUBROUTINE BasisTOGrid_Basis_rapide


  SUBROUTINE BasisTOGrid_Basis_rapide1(G,B,Basis)
  USE UtilLib_m
  USE NDindex_m
    !Logical,           parameter          :: debug = .true.
    Logical,          parameter             :: debug = .false.
    TYPE(Basis_t),  intent(in)              :: Basis
    Real (kind=Rk),  intent(in),target      :: B(:) !Vector on base,
    Real (kind=Rk),  intent(inout),target   :: G(:) !Vector on the grid, out
    Real (kind=Rk), pointer                 :: BBG(:,:,:)
    Real (kind=Rk), pointer                 :: BBB(:,:,:)
    !Real (kind=Rk) ,allocatable,target      :: GB1(:),GB2(:)
    integer , allocatable                   :: Ib3(:),Iq1(:),Iq2(:),ib1(:),ib2(:),iq3(:)
    Real (kind=Rk), pointer                 :: GBB(:,:,:)
    Real (kind=Rk) ,allocatable,target      :: GBB1(:),GGB2(:)

    Integer                                 :: ib,iq,nq,nb,inb,Ndim
    Integer                                 :: jb,jb1,jb2

    IF(debug)THEN
     write(out_unitp,*) 'BEGINNING BasisTOGrid_Basis'
     write(out_unitp,*) 'intent(in) :: B(:)',B
     Call Write_Basis(Basis)
     flush(out_unitp)
    END IF

    IF(.NOT. Basis_IS_Allocated(Basis))THEN
     write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
     write(out_unitp,*) " the basis is not Allocated."
     STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
    END IF

    IF(size(B) /= Basis%nb)THEN
     write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
     write(out_unitp,*) ' the size of B is different from nb.'
     write(out_unitp,*) ' size(B), Basis%nb',size(B),Basis%nb
     STOP 'ERROR in BasisTOGrid_Basis: wrong B size.'
    END IF

    IF(size(G) /= Basis%nq)THEN
     write(out_unitp,*) ' ERROR in GridTOBasis_Basis'
     write(out_unitp,*) ' the size of G is different from nq.'
     write(out_unitp,*) ' size(G), Basis%nq',size(G),Basis%nq
     STOP 'ERROR in BasisTOGrid_Basis: wrong G size..'
    END IF


    IF(Allocated(Basis%tab_basis)) THEN
      Ndim = size(Basis%tab_basis)

      Call Calc_indice0( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)

      Allocate(GBB1(iq1(1)*iq2(1)*ib3(1)))

      BBB(1:iq1(1),1:ib2(1),1:ib3(1)) => B
      GBB(1:iq1(1),1:iq2(1),1:ib3(1)) => GBB1

      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(1))

      DO inb = 2,Ndim-1

        Allocate(GGB2(iq1(inb)*iq2(inb)*ib3(inb)))

        BBB(1:iq1(Inb),1:ib2(inb),1:ib3(inb)) => GBB1

        GBB(1:iq1(inb),1:iq2(inb),1:ib3(inb)) => GGB2

        Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(inb))

        GBB1=GGB2

        Deallocate(GGB2)

      END DO

      BBB(1:iq1(Ndim),1:ib2(Ndim),1:ib3(Ndim)) => GBB1

      GBB(1:iq1(Ndim),1:iq2(Ndim),1:ib3(Ndim)) => G

      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(Ndim))
      Deallocate(Ib1,Iq1,Iq2,Ib2,Ib3,Iq3)


    ELSE
      DO iq = 1,Basis%nq
       G(iq) = ZERO
       DO ib = 1,Basis%nb
        G(iq) = G(iq) + Basis%d0gb(iq,ib) * B(ib)
       END DO
      END DO
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'intent(OUTIN) :: G(:)',G
      write(out_unitp,*) 'END BasisTOGrid_Basis_rapide'
      flush(out_unitp)
    END IF
  END SUBROUTINE BasisTOGrid_Basis_rapide1


  SUBROUTINE GridTOBasis_Basis(B,G,Basis)
  USE UtilLib_m
    !Logical,          parameter    :: debug = .true.
    Logical,         parameter     :: debug = .false.
    TYPE(Basis_t),    intent(in)    :: Basis
    Real(kind=Rk),    intent(in)    :: G(:) !Vector on the grid, at the entrance
    Real(kind=Rk),    intent(inout) :: B(:) !Vector on base, out
    Logical                         :: Endloop_q
    Logical                         :: Endloop_b
    Real(kind=Rk)                   :: WT,W
    Integer,        allocatable     :: tab_iq(:)
    Integer,        allocatable     :: tab_ib(:)
    Integer                         :: ib,iq,iq1,iq2,inb
    Integer                         :: jb,ib1,ib2,jb1,jb2

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING GridTOBasis_Basis'
      write(out_unitp,*) 'intent(in) :: G(:)',G
      !Call Write_Basis(Basis)
      flush(out_unitp)
    END IF

    IF (.NOT. Basis_IS_Allocated(Basis)) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
      write(out_unitp,*) " the basis is not Allocated."
      STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
    END IF

    IF (size(B) /= Basis%nb) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
      write(out_unitp,*) ' the size of G is different from nb.'
      write(out_unitp,*) ' size(B), Basis%nb',size(B),Basis%nb
      STOP 'ERROR in GridTOBasis_Basis: wrong B size.'
    END IF

    IF (size(G) /= Basis%nq) THEN
      write(out_unitp,*) ' ERROR in GridTOBasis_Basis'
      write(out_unitp,*) ' the size of G is different from nq.'
      write(out_unitp,*) ' size(G), Basis%nq',size(G),Basis%nq
      STOP 'ERROR in GridTOBasis_Basis: wrong G size'
    END IF

    IF(Allocated(Basis%tab_basis)) THEN
     Allocate(Tab_ib(size(Basis%tab_basis)))
     Allocate(Tab_iq(size(Basis%tab_basis)))
     Call Init_tab_ind(Tab_ib,Basis%NDindexb)
     Ib = 0
     DO
      Ib = Ib+1
      Call increase_NDindex(Tab_ib,Basis%NDindexb,Endloop_b)
      IF (Endloop_b) exit
      B(ib) = ZERO
      Call Init_tab_ind(Tab_iq,Basis%NDindexq)
      Iq = 0
      DO
        Iq = Iq+1
        Call increase_NDindex(Tab_iq,Basis%NDindexq,Endloop_q)
        IF (Endloop_q) exit
        WT = 1
        W  = 1
        DO inb = 1,size(Basis%tab_basis)
         WT = WT * Basis%tab_basis(inb)%w(tab_iq(inb))
         W  = W  * Basis%tab_basis(inb)%d0gb(tab_iq(inb),tab_ib(inb))
        END DO
        B(ib) = B(ib) + W * WT * G(iq)
      END DO
    END DO
     DeAllocate(Tab_ib)
     DeAllocate(Tab_iq)
    ELSE
     DO ib = 1,Basis%nb
       B(ib) = ZERO
       DO iq = 1,Basis%nq
         B(ib) = B(ib) + Basis%d0gb(iq,ib) * Basis%w(iq) * G(iq)
       END DO
     END DO
    END IF

    IF (debug) THEN
     ! write(out_unitp,*) 'intent(OUTIN) :: B(:)',B
     write(out_unitp,*) 'END GridTOBasis_Basis'
     flush(out_unitp)
    END IF
  END SUBROUTINE GridTOBasis_Basis



  SUBROUTINE GridTOBasis_Basis_rapide(B,G,Basis)
  USE UtilLib_m
    !Logical,          parameter    :: debug = .true.
    Logical,         parameter     :: debug = .false.
    TYPE(Basis_t),   intent(in),target        :: Basis
    Real (kind=Rk),  intent(in) ,target       :: G(:)
    Real (kind=Rk),  intent(inout),target     :: B(:)
    Real (kind=Rk),  pointer                  :: BBB(:,:,:),BGG(:,:,:)
    Real(kind=Rk) ,  allocatable  ,target     :: BGG1(:),BGG2(:)
    Real (kind=Rk),  pointer                  :: GGG(:,:,:),BBG(:,:,:)
    Integer                                   :: ib,i1,i3,inb,Ndim,iq
    Integer,         allocatable              :: Ib1(:),Ib2(:),Iq3(:),Iq1(:),Iq2(:),Ib3(:)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING GridTOBasis_Basis_rapide'
      write(out_unitp,*) 'intent(in) :: G(:)',G
      !Call Write_Basis(Basis)
      flush(out_unitp)
    END IF

    IF (.NOT. Basis_IS_Allocated(Basis)) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basi_rapides'
      write(out_unitp,*) " the basis is not Allocated."
      STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
    END IF

    IF (size(B) /= Basis%nb) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basis_rapide'
      write(out_unitp,*) ' the size of G is different from nb.'
      write(out_unitp,*) ' size(B), Basis%nb',size(B),Basis%nb
      STOP 'ERROR in GridTOBasis_Basis: wrong B size.'
    END IF

    IF (size(G) /= Basis%nq) THEN
      write(out_unitp,*) ' ERROR in GridTOBasis_Basis'
      write(out_unitp,*) ' the size of G is different from nq.'
      write(out_unitp,*) ' size(G), Basis%nq',size(G),Basis%nq
      STOP 'ERROR in GridTOBasis_Basis: wrong G size'
    END IF

    IF(Allocated(Basis%tab_basis)) THEN

      Allocate(BGG1(Basis%tab_basis(1)%nb*Basis%tab_basis(2)%nq*Basis%tab_basis(3)%nq))

      BGG1(:) = ZERO

      GGG(1:1,1:Basis%tab_basis(1)%nq,1:Basis%tab_basis(2)%nq*Basis%tab_basis(3)%nq)   => G
      BGG(1:1,1:Basis%tab_basis(1)%nb,1:Basis%tab_basis(2)%nq*Basis%tab_basis(3)%nq)   => BGG1


      Call GridTOBasis_1D(BGG,GGG,Basis%tab_basis(1))

      Allocate(BGG2(Basis%tab_basis(1)%nb*Basis%tab_basis(2)%nb*Basis%tab_basis(3)%nq))

      BGG2(:) = ZERO

      GGG( 1:Basis%tab_basis(1)%nb,1:Basis%tab_basis(2)%nq,1:Basis%tab_basis(3)%nq)    => BGG1
      BGG( 1:Basis%tab_basis(1)%nb,1:Basis%tab_basis(2)%nb,1:Basis%tab_basis(3)%nq)    => BGG2

      Call GridTOBasis_1D(BGG,GGG,Basis%tab_basis(2))

      B(:) = ZERO

      GGG(1:Basis%tab_basis(1)%nb*Basis%tab_basis(2)%nb,1:Basis%tab_basis(3)%nq,1:1)   => BGG2
      BGG(1:Basis%tab_basis(1)%nb*Basis%tab_basis(2)%nb,1:Basis%tab_basis(3)%nb,1:1)   => B

      Call GridTOBasis_1D(BGG,GGG,Basis%tab_basis(3))


    ELSE
      DO ib = 1,Basis%nb
       B(ib) = ZERO
       DO iq = 1,Basis%nq
         B(ib) = B(ib) + Basis%d0gb(iq,ib) * Basis%w(iq) * G(iq)
       END DO
      END DO
    END IF

    IF (debug) THEN
     write(out_unitp,*) 'END GridTOBasis_Basis_rapide'
     flush(out_unitp)
    END IF
  END SUBROUTINE GridTOBasis_Basis_rapide

  SUBROUTINE GridTOBasis_Basis_rapide1(B,G,Basis)
  USE UtilLib_m
    !Logical,          parameter    :: debug = .true.
    Logical,         parameter     :: debug = .false.
    TYPE(Basis_t),   intent(in),target        :: Basis
    Real (kind=Rk),  intent(in) ,target       :: G(:)
    Real (kind=Rk),  intent(inout),target     :: B(:)
    Real (kind=Rk),  pointer                  :: BBB(:,:,:),GGB(:,:,:)
    Real(kind=Rk) ,  allocatable  ,target     :: BGG1(:),BGG2(:)
    Real (kind=Rk),  pointer                  :: GGG(:,:,:)!,BBG(:,:,:)
    Integer                                   :: ib,i1,i3,inb,Ndim,iq
    Integer,         allocatable              :: Ib1(:),Ib2(:),Iq3(:),Iq1(:),Iq2(:),Ib3(:)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING GridTOBasis_Basis_rapide'
      write(out_unitp,*) 'intent(in) :: G(:)',G
      !Call Write_Basis(Basis)
      flush(out_unitp)
    END IF

    IF (.NOT. Basis_IS_Allocated(Basis)) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basi_rapides'
      write(out_unitp,*) " the basis is not Allocated."
      STOP "ERROR BasisTOGrid_Basis: the basis is not Allocated."
    END IF

    IF (size(B) /= Basis%nb) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basis_rapide'
      write(out_unitp,*) ' the size of G is different from nb.'
      write(out_unitp,*) ' size(B), Basis%nb',size(B),Basis%nb
      STOP 'ERROR in GridTOBasis_Basis: wrong B size.'
    END IF

    IF (size(G) /= Basis%nq) THEN
      write(out_unitp,*) ' ERROR in GridTOBasis_Basis'
      write(out_unitp,*) ' the size of G is different from nq.'
      write(out_unitp,*) ' size(G), Basis%nq',size(G),Basis%nq
      STOP 'ERROR in GridTOBasis_Basis: wrong G size'
    END IF

    IF(Allocated(Basis%tab_basis)) THEN

      Ndim = size(Basis%tab_basis)
      Call Calc_indice0( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
      Allocate(BGG1(Ib1(1)*Ib2(1)*Iq3(1)))

      BGG1=ZERO

      GGG(1:Ib1(1),1:Iq2(1),1:Iq3(1))   => G
      GGB(1:Ib1(1),1:Ib2(1),1:Iq3(1))    => BGG1

      Call GridTOBasis_1D(GGB,GGG,Basis%tab_basis(1))

      DO inb = 2,Ndim-1

        Allocate(BGG2(Ib1(inb)*Ib2(inb)*Iq3(inb)))

        BGG2=ZERO

        GGG( 1:Ib1(inb),1:Iq2(inb),1:Iq3(inb))    => BGG1
        GGB( 1:Ib1(inb),1:Ib2(inb),1:Iq3(inb))    => BGG2

        Call GridTOBasis_1D(GGB,GGG,Basis%tab_basis(inb))

        BGG1=BGG2

        deallocate(BGG2)
      END DO

      B(:) = ZERO

      GGG(1:Ib1(Ndim),1:Iq2(Ndim),1:Iq3(Ndim)) => BGG1
      GGB(1:Ib1(Ndim),1:Ib2(Ndim),1:Iq3(Ndim)) => B

      Call GridTOBasis_1D(GGB,GGG,Basis%tab_basis(Ndim))

      Deallocate (Iq1,Iq2,Iq3,Ib1,Ib2,Ib3)
    ELSE
      DO ib = 1,Basis%nb
       B(ib) = ZERO
       DO iq = 1,Basis%nq
         B(ib) = B(ib) + Basis%d0gb(iq,ib) * Basis%w(iq) * G(iq)
       END DO
      END DO
    END IF

    IF (debug) THEN
     write(out_unitp,*) 'END GridTOBasis_Basis_rapide'
     flush(out_unitp)
    END IF
  END SUBROUTINE GridTOBasis_Basis_rapide1


  SUBROUTINE Calc_dngg_grid(Basis)
  USE UtilLib_m
    TYPE(Basis_t), intent(inout)    :: Basis
    Real(kind=Rk), allocatable      :: d0bgw(:,:)
    Integer                         :: ib,i1,i2
    !Logical,          parameter    :: debug = .true.
    Logical,         parameter    ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Calc_dngg_grid'
      Call Write_Basis(Basis)
      flush(out_unitp)
    END IF

    Allocate(Basis%d1gg(Basis%nq,Basis%nq,Basis%ndim))
    Allocate(Basis%d2gg(Basis%nq,Basis%nq,Basis%ndim,Basis%ndim))

    d0bgw = transpose(Basis%d0gb)
    DO ib=1,Basis%nb
       d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
    END DO

    IF (debug) THEN
      Call Write_RMat(d0bgw(:,:),out_unitp,5,name_info='d0bgw')
      write(out_unitp,*)
    END IF

    DO i1=1,Basis%ndim
        Basis%d1gg(:,:,i1)   = matmul(Basis%d1gb(:,:,i1),d0bgw)

       DO i2=1,Basis%ndim
         Basis%d2gg(:,:,i1,i2) = matmul(Basis%d2gb(:,:,i1,i2),d0bgw)
       END DO
    END DO
!STOP 'FFFFFF'
    !Call Write_RMat(Basis%d1gg(:,:,1),out_unitp,5,name_info='d1gg')
    !write(out_unitp,*)
    !Call Write_RMat(Basis%d2gg(:,:,1,1),out_unitp,5,name_info='d2gg')

    IF (debug) THEN
      Call Write_Basis(Basis)
      write(out_unitp,*) 'END Calc_dngg_grid'
      flush(out_unitp)
    END IF
    deAllocate(d0bgw)
  END SUBROUTINE Calc_dngg_grid

  SUBROUTINE Calc_dngg_grid_old(Basis)
  USE UtilLib_m
    TYPE(Basis_t), intent(inout)    :: Basis
    Real(kind=Rk), allocatable      :: d0bgw(:,:)
    Integer                         :: ib
    !Logical,          parameter    :: debug = .true.
    Logical,         parameter    ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Calc_dngg_grid_old'
      Call Write_Basis(Basis)
      flush(out_unitp)
    END IF

    Allocate(Basis%d1gg(Basis%nq,Basis%nq,1))
    Allocate(Basis%d2gg(Basis%nq,Basis%nq,1,1))

    d0bgw = transpose(Basis%d0gb)
    DO ib=1,Basis%nb
       d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
    END DO

    IF (debug) THEN
      Call Write_RMat(d0bgw(:,:),out_unitp,5,name_info='d0bgw')
      write(out_unitp,*)
    END IF

    Basis%d1gg(:,:,1)   = matmul(Basis%d1gb(:,:,1),d0bgw)
    Basis%d2gg(:,:,1,1) = matmul(Basis%d2gb(:,:,1,1),d0bgw)

    !Call Write_RMat(Basis%d1gg(:,:,1),out_unitp,5,name_info='d1gg')
    !write(out_unitp,*)
    !Call Write_RMat(Basis%d2gg(:,:,1,1),out_unitp,5,name_info='d2gg')

    IF (debug) THEN
      Call Write_Basis(Basis)
      write(out_unitp,*) 'END Calc_dngg_grid_old'
      flush(out_unitp)
    END IF
    deAllocate(d0bgw)
  END SUBROUTINE Calc_dngg_grid_old

  SUBROUTINE Test_Passage(Basis)
  USE UtilLib_m
    TYPE(Basis_t),    intent(in)    :: Basis
  !  Logical,          parameter    :: debug = .true.
    Logical,         parameter    ::debug = .false.
    Real(kind=Rk),    allocatable   :: G1(:),B1(:)
    Real(kind=Rk),    allocatable   :: G2(:),B2(:)
    Real(kind=Rk),    allocatable   :: B(:),Delta_g(:),Delta_b(:)
    Real(kind=Rk),    parameter    :: eps = TEN**(-TEN)
    Integer                         :: iq

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Test_Passage'
      Call Write_Basis(Basis)
      flush(out_unitp)
    END IF

    Allocate(B(Basis%nb))
    Allocate(Delta_g(Basis%nq))
    Allocate(G1(Basis%nq))
    Allocate(G2(Basis%nq))

    G1(:)=ONE
    Call GridTOBasis_Basis(B,G1,Basis)
    Call BasisTOGrid_Basis(G2,B,Basis)
    Delta_g(:) = ABS(G1(:) - G2(:))

    IF (eps.gt.maxval(Delta_g(:)) ) THEN
      Write(out_unitp,*) ' Best translate for Test_Grid_Basis_Grid '
      Write(out_unitp,*) maxval(Delta_g(:))
    ELSE
      Write(out_unitp,*) maxval(Delta_g(:))
      STOP 'Bad translate for Test_Grid_Basis_Grid'
    END IF

    IF (Basis%nq>=Basis%nb ) THEN
      Allocate(Delta_b(Basis%nb))
      Allocate(B1(Basis%nb))
      Allocate(B2(Basis%nb))
      B1(:)=ONE
      Call BasisTOGrid_Basis(G1,B1,Basis)
      Call GridTOBasis_Basis(B2,G1,Basis)
      Delta_b(:) = ABS(B1(:) - B2(:))
      IF (eps.gt.maxval(Delta_b(:)) ) THEN
        Write(out_unitp,*) ' Best translate for Basis_Grid_Basis '
        Write(out_unitp,*) maxval(Delta_b(:))
      ELSE
        STOP 'Bad translate for Basis_Grid_Basis'
      END IF
    END IF
    IF (debug) THEN
      write(out_unitp,*) 'END Test_Passage'
      Call Write_Basis(Basis)
      flush(out_unitp)
    END IF
END SUBROUTINE Test_Passage

SUBROUTINE Scale_Basis(Basis)
USE UtilLib_m
    TYPE(Basis_t),       intent(inout)  :: Basis
    Real(kind=Rk)                      :: x0,sx

    IF(Basis%Ndim==1)THEN
      x0 = Basis%Q0
      sx = Basis%scaleQ
    ELSE
      x0 = Zero
      sx = One
    END IF
    IF (abs(sx) > ONETENTH**6 .AND. Basis_IS_Allocated(Basis)) THEN

      Basis%x(:,1) = x0 + Basis%x(:,1) / sx
      Basis%w(:) =      Basis%w(:) / sx

      Basis%d0gb(:,:)     = Basis%d0gb(:,:)     * sqrt(sx)
      Basis%d1gb(:,:,:)   = Basis%d1gb(:,:,:)   * sqrt(sx)*sx
      Basis%d2gb(:,:,:,:) = Basis%d2gb(:,:,:,:) * sqrt(sx)*sx*sx
    ELSE
      write(out_unitp,*) ' ERROR in Scale_Basis'
      write(out_unitp,*) ' sx is too small  or ...'
      write(out_unitp,*) ' the basis is not Allocated.'
      write(out_unitp,*) 'Sx', Sx
      write(out_unitp,*) 'alloc', Basis_IS_Allocated(Basis)
      STOP 'ERROR in Scale_Basis'
    END IF
  END SUBROUTINE Scale_Basis

END MODULE Basis_m
