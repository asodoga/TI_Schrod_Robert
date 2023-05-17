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
  USE NDindex_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Basis_t,Basis_IS_Allocated,BasisTOGrid_Basis,GridTOBasis_Basis
  PUBLIC :: Test_Passage,Calc_dngg_grid,Basis_IS_Allocatedtot,write_basis
  PUBLIC :: BasisTOGrid_Basis_rapide,GridTOBasis_Basis_rapide,Read_Construct_Basis
  PUBLIC :: BasisTOGrid_Basis_rapide1,GridTOBasis_Basis_rapide1,NDSmolyak_t
  PUBLIC :: Read_Basis_old,Calc_n_smol,IND_Smolyak1_t,Calc_n_smol_1


  TYPE :: IND_Smolyak1_t
    Integer, allocatable              :: Nb_smol_1(:)
  END TYPE  IND_Smolyak1_t

  TYPE :: NDSmolyak_t
    Integer                       :: NSB   = 0
    Integer                       :: NSQ   = 0
    Integer                       :: Ndim  = 0
    Integer                       :: Maxtermsmol  = 0
    Integer, allocatable          :: Li_smo(:,:)
    Integer, allocatable          :: Nqli(:,:)
    Integer, allocatable          :: Nbli(:,:)
    Integer, allocatable          :: DSmol(:)
    Integer, allocatable          :: Nblp(:)
    Integer, allocatable          :: Nqlp(:)
    !Integer, allocatable          :: Tags1(:)
    Integer, allocatable          :: Sum_Nqlp(:)
    Integer, allocatable          :: Sum_Nblp(:)
    Integer, allocatable          :: NDend_smol_q(:,:)
    Integer, allocatable          :: NDend_smol_b(:,:)
    TYPE (IND_Smolyak1_t),allocatable :: SMOL_1(:)
    TYPE (IND_Smolyak1_t),allocatable :: SMOL_3(:)
  END TYPE NDSmolyak_t

  TYPE :: Basis_t
    Integer                      :: nb_basis   = 0
    Integer                      :: nb         = 0
    Integer                      :: nq         = 0
    Integer, allocatable         :: A_smol(:),B_smol(:)
    Integer, allocatable         :: Tags1(:)
    Integer, allocatable         :: Tags2(:)
    Integer, allocatable         :: Ind_map(:)
    Real(kind=Rk)                :: A,B,scaleQ,Q0
    Character(len=:),allocatable :: Basis_name
    Real(kind=Rk),   allocatable :: x(:)
    Real(kind=Rk),   allocatable :: w(:)
    Real(kind=Rk),   allocatable :: d0gb(:,:)      ! basis s d0gb(nq,nb)
    Real(kind=Rk),   allocatable :: d1gb(:,:,:)    ! basis s d2gb(nq,nb,1)
    Real(kind=Rk),   allocatable :: d1gg(:,:,:)    ! basis s d2gg(nq,nq,1)
    Real(kind=Rk),   allocatable :: d2gb(:,:,:,:)  ! basis s d2gb(nq,nb,1,1)
    Real(kind=Rk),   allocatable :: d2gg(:,:,:,:)  ! basis s d2gg(nq,nq,1,1)
    TYPE(NDindex_t), allocatable :: NDindexq_smol(:)
    TYPE(NDindex_t), allocatable :: NDindexb_smol(:)
    TYPE(NDindex_t)              :: NDindexq
    TYPE(NDindex_t)              :: NDindexb
    TYPE(NDindex_t)              :: NDindexl
    TYPE(NDindex_t)              :: NDindex_smol_q
    TYPE(NDindex_t)              :: NDindex_smol_b
    TYPE(NDSmolyak_t)            :: Smolyak
    TYPE(NDSmolyak_t),allocatable:: Smolyak2(:)
    TYPE(NDSmolyak_t),allocatable:: Smolyak1(:)
    TYPE (Basis_t),  allocatable :: tab_basis(:)
    TYPE (Basis_t),  allocatable :: tab_Smolyak(:)
    TYPE (Basis_t),  allocatable :: tab_in_Smolyak(:)
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

      write(out_unitp,*)"Aet B",Basis%A,Basis%B
      write(out_unitp,*) 'Scaleq',Basis%Scaleq
      write(out_unitp,*) 'Q0' ,Basis%Q0
      IF (.NOT.Allocated(Basis%x)) THEN
       write(out_unitp,*)' Basis table x is not Allocated.'
      ELSE
        Call Write_RVec(Basis%x,out_unitp,5,name_info='x')
      END IF
      write(out_unitp,*)
      IF (.NOT.Allocated(Basis%W)) THEN
        write(out_unitp,*)' Basis table w is not Allocated.'
      ELSE
        Call Write_RVec(Basis%w,out_unitp,5,name_info='w')
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
    Logical,        parameter               :: debug = .true.
    !Logical,       parameter               :: debug = .false.
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

  SUBROUTINE Construct_basis_new (Basis)
   USE UtilLib_m
     Logical,              parameter      :: debug = .true.
     !Logical,             parameter       ::debug = .false.
     TYPE(Basis_t),  intent(inout)        :: Basis

     Call string_uppercase_TO_lowercase(Basis%Basis_name)
     SELECT CASE (Basis%Basis_name)
     CASE ('boxab')
     Call Construct_Basis_Sin(Basis)
     Basis%Q0      = Basis%A
     Basis%scaleQ  = pi/(Basis%B-Basis%A)
     CASE ('herm','ho')
     Call Construct_Basis_Ho(Basis)
     CASE default
     STOP 'ERROR in Read_Basis: no default basis.'
     END SELECT

     Call Scale_Basis(Basis)
     Call Calc_dngg_grid(Basis)
     Call CheckOrtho_Basis(Basis,nderiv=2)

  END SUBROUTINE Construct_basis_new


  FUNCTION Calc_n_smol(A,B,l)
   Integer, intent(in)     :: A
   Integer, intent(in)     :: B
   Integer, intent(in)     :: L
   Integer                 :: Calc_n_smol

    Calc_n_smol = A + B*L

  END FUNCTION Calc_n_smol


  SUBROUTINE Calc_n_smol_1(Basis,A,B,J)

   TYPE(Basis_t),  intent(inout):: Basis
   Integer                      :: I,s,n
   Integer, intent(in)          :: A
   Integer, intent(in)          :: B
   Integer, intent(in)          :: J

   !DO J=1, Basis%Smolyak%Ndim
   Allocate(Basis%Smolyak1(J)%Smol_1(Basis%NDindexl%L+1))
   DO I=1,Basis%NDindexl%L+1
    !write(out_unitp,*)
     IF (I==1) THEN
       n = Calc_n_smol(A,B,I-1)
       Allocate(Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(n))
       DO s=1,n
         Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(s)=s

         write(out_unitp,*) ' Nb_smol_1(s)', Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(s)
       END DO
     ELSE
       n = Calc_n_smol(A,B,I-1)-Calc_n_smol(A,B,I-2)
       Allocate(Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(n))
       Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(1) = Calc_n_smol(A,B,I-2) +1
       write(out_unitp,*) ' Nb_smol_1(s)',Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(1)
       !write(out_unitp,*)
       DO s = 2, n
        Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(s) = Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(s-1)+1
        write(out_unitp,*) ' Nb_smol_1(s)',Basis%Smolyak1(J)%Smol_1(I)%Nb_smol_1(s)
       END DO
     END IF
     write(out_unitp,*)
   END DO
!write(out_unitp,*)
  END SUBROUTINE Calc_n_smol_1



  SUBROUTINE Mapping_2 (Basis )
    TYPE(Basis_t),intent(inout)  :: Basis
    !Integer, intent(in)          :: Tab_ind1(:)
    !Integer, intent(in)          :: J
    !Integer, intent(inout)       :: Tags2(:,:)
    Integer                      :: I_smol,n,i,I_smol1
    Integer, allocatable         :: Tab_ind(:)
    Integer, allocatable         :: Tab_li(:)
    Integer, allocatable         :: Tags2(:,:)
    Integer, allocatable         :: Ndend(:,:)
    logical                      :: Endloop

    Allocate(Tab_li(Basis%Smolyak%Ndim))

    Allocate(Ndend(Basis%Smolyak%Maxtermsmol,Basis%Smolyak%Ndim))
    Call Init_tab_ind(Tab_li,Basis%NDindexl)
    I_smol = 0
    DO
      I_smol = I_smol+1
      CALL increase_NDindex(Tab_li,Basis%NDindexl,Endloop)
      Allocate(Basis%tab_Smolyak(I_smol)%Tags1(Basis%Smolyak%Ndim+1))
      DO i=1,Basis%Smolyak%Ndim
        Basis%tab_Smolyak(I_smol)%Tags1(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))
      END DO
      Basis%tab_Smolyak(I_smol)%Tags1(Basis%Smolyak%Ndim+1) = I_smol

      Write(*,*) 'Tab',Basis%tab_Smolyak(I_smol)%Tags1

      DO i=1,Basis%Smolyak%Ndim
        Ndend(I_smol,i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))
      END DO

      IF (Endloop) exit
    END DO

    Allocate(Tab_ind(Basis%Smolyak%Ndim))

    DO I_smol=1,Basis%Smolyak%Maxtermsmol
      Allocate(Tags2(Basis%tab_Smolyak(I_smol)%nb,Basis%Smolyak%Ndim+1))
      Tab_ind(:) = 1
      Tab_ind(1) = 0
      Do n=1,Basis%tab_Smolyak(I_smol)%nb
        Tab_ind(1) = Tab_ind(1)+1
        DO i=1,Basis%Smolyak%Ndim-1
          IF(Tab_ind(i)  > Ndend(I_smol,i)) THEN
           Tab_ind(i+1) = Tab_ind(i+1)+1
           Tab_ind(i)   = 1
          END IF
        END DO


        DO I_smol1=1,I_smol
           IF( ALL(Tab_ind(1:Basis%Smolyak%Ndim)==Basis%tab_Smolyak(I_smol1)%Tags1(1:Basis%Smolyak%Ndim)))THEN

             Tags2(n,:)=Basis%tab_Smolyak(I_smol1)%Tags1(:)

           END IF
        END DO

        Write(*,*) 'mapping',Tags2(n,:)

      END DO
      Deallocate(Tags2)
      Write(*,*) '-------------------------------------------------------'
    END DO

  END SUBROUTINE Mapping_2


  SUBROUTINE Mapping_2G (Basis )
    Logical,             parameter      :: debug = .true.
    !Logical,             parameter       :: debug = .false.
    TYPE(Basis_t),intent(inout)  :: Basis
    Integer                      :: I_smol,n,i,I_smol1,compt
    Integer, allocatable         :: Tab_ind(:)
    Integer, allocatable         :: Tab_li(:)
    Integer, allocatable         :: Tags2(:,:)
    Integer, allocatable         :: Ndend(:,:)
    Integer, allocatable         :: Ndbeging(:)
    Integer                      :: Nbinter,I_smol2,I_compt
    Integer, allocatable         :: Nbinter1(:)
    Integer, allocatable         :: Tags(:)
    logical                      :: Endloop


    Allocate(Tab_li(Basis%Smolyak%Ndim))
    Allocate(Ndbeging(Basis%Smolyak%Ndim))
    Allocate(Nbinter1(Basis%Smolyak%Ndim))
    Allocate(Ndend(Basis%Smolyak%Maxtermsmol,Basis%Smolyak%Ndim))

    Call Init_tab_ind(Tab_li,Basis%NDindexl)

    I_smol = 0
    Compt  = 0

    DO
      I_smol = I_smol+1

      CALL increase_NDindex(Tab_li,Basis%NDindexl,Endloop)

      Allocate(Basis%tab_Smolyak(I_smol)%Tags1(Basis%Smolyak%Ndim))

      DO i=1,Basis%Smolyak%Ndim
        Ndend(I_smol,i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))
      END DO

      DO i = 1,Basis%Smolyak%Ndim
        IF(Tab_li(i)==0)THEN
         Ndbeging(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))
         Nbinter1(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))
        ELSE
         Ndbeging(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i)-1) +1
         Nbinter1(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i)) &
                -Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i)-1)
        END IF
      END DO

      Nbinter  = Product(Nbinter1(:))

      Basis%tab_Smolyak(I_smol)%Tags1(1:Basis%Smolyak%Ndim) = Ndbeging(:)
      Basis%tab_Smolyak(I_smol)%Tags1(1) = Ndbeging(1)-1

      DO  n = 1,Nbinter
        Basis%tab_Smolyak(I_smol)%Tags1(1) = Basis%tab_Smolyak(I_smol)%Tags1(1)+1
        DO i = 1,Basis%Smolyak%Ndim-1
          IF(Basis%tab_Smolyak(I_smol)%Tags1(i)  >  Ndend(I_smol,i)) THEN
           Basis%tab_Smolyak(I_smol)%Tags1(i+1) = Basis%tab_Smolyak(I_smol)%Tags1(i+1)+1
           Basis%tab_Smolyak(I_smol)%Tags1(i)   = Ndbeging(i)
          END IF
        END DO
        Compt=Compt+1

        Write(*,*) 'Tab',Basis%tab_Smolyak(I_smol)%Tags1
      END DO
      Deallocate(Basis%tab_Smolyak(I_smol)%Tags1)

      IF (Endloop) exit

    END DO

    Allocate(Basis%tab_in_Smolyak(Compt))


    Call Init_tab_ind(Tab_li,Basis%NDindexl)

    I_smol2 = 0
    Compt   = 0

    DO
      I_smol2 = I_smol2 + 1

      CALL increase_NDindex(Tab_li,Basis%NDindexl,Endloop)

      Allocate(Tags(Basis%Smolyak%Ndim+1))


      DO i = 1,Basis%Smolyak%Ndim
        IF(Tab_li(i)==0)THEN
         Ndbeging(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))
         Nbinter1(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))
        ELSE
         Ndbeging(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i)-1) +1
         Nbinter1(i) = Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i))&
                -Calc_n_smol(Basis%A_smol(i),Basis%B_smol(i),Tab_li(i)-1)
        END IF
      END DO

      Nbinter  = Product(Nbinter1(:))

      Tags(1:Basis%Smolyak%Ndim) = Ndbeging(:)
      Tags(1) = Ndbeging(1)-1

      DO  n = 1,Nbinter
        Tags(1) = Tags(1)+1
        DO i = 1,Basis%Smolyak%Ndim-1
          IF(Tags(i)  >  Ndend(I_smol2,i)) THEN
           Tags(i+1) =  Tags(i+1)+1
           Tags(i)   =  Ndbeging(i)
          END IF
        END DO

        Compt=Compt+1

        Allocate(Basis%tab_in_Smolyak(Compt)%Tags1(Basis%Smolyak%Ndim+1))

        Tags(Basis%Smolyak%Ndim+1)        = Compt
        Basis%tab_in_Smolyak(Compt)%Tags1 = Tags

      END DO

      Deallocate(Tags)

      IF (Endloop) exit
    END DO

    Allocate(Tab_ind(Basis%Smolyak%Ndim))

    DO I_smol=1,Basis%Smolyak%Maxtermsmol
      Allocate(Tags2(Basis%tab_Smolyak(I_smol)%nb,Basis%Smolyak%Ndim+1))
      Allocate(Basis%tab_Smolyak(I_smol)%Ind_map(Basis%tab_Smolyak(I_smol)%nb))
      Tab_ind(:) = 1
      Tab_ind(1) = 0
      Do n=1,Basis%tab_Smolyak(I_smol)%nb
        Tab_ind(1) = Tab_ind(1)+1
        DO i=1,Basis%Smolyak%Ndim-1
          IF(Tab_ind(i)  > Ndend(I_smol,i)) THEN
           Tab_ind(i+1) = Tab_ind(i+1)+1
           Tab_ind(i)   = 1
          END IF
        END DO

        DO I_compt = 1,Compt

           IF( ALL(Tab_ind(1:Basis%Smolyak%Ndim) == Basis%tab_in_Smolyak(I_compt)%Tags1(1:Basis%Smolyak%Ndim)))THEN
             Tags2(n,:) = Basis%tab_in_Smolyak(I_compt)%Tags1(:)
             Basis%tab_Smolyak(I_smol)%Ind_map(n)= I_compt
             exit
           END IF

        END DO
        Write(*,*) 'Ind_map',Basis%tab_Smolyak(I_smol)%Ind_map(n)
        !Write(*,*) 'mapping',Tags2(n,:)
      END DO
      Deallocate(Tags2)
      Write(*,*) '-------------------------------------------------------'
    END DO

  END SUBROUTINE Mapping_2G

  SUBROUTINE Read_Construct_Basis(Basis,nio)
   USE UtilLib_m
    Logical,             parameter      :: debug = .true.
   !Logical,             parameter       :: debug = .false.
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
    Integer                             :: err_io,nb,nq,i,j,nb_basis,LB,NDim,inb,D
    Integer                             :: A_smol,B_smol,S,n
    Character (len=Name_len)            :: name
    Real(kind=Rk)                       :: A,B,scaleQ,Q0,d0,d2,X1,W1

    Call Read_Basis_simple(name,nb_basis,LB,A_smol,B_smol,nb,nq,A,B,scaleQ,Q0,nio)
    Basis%Basis_name        = name

    IF (nb_basis > 1) THEN

      NDim = nb_basis
      Basis%NDindexq%sys_type = name
      Basis%NDindexb%sys_type = name
      Basis%NDindexl%sys_type = name
      Basis%NDindexl%L        = LB

      Call string_uppercase_TO_lowercase(Basis%Basis_name)
      Call string_uppercase_TO_lowercase(Basis%NDindexl%sys_type)
      Call string_uppercase_TO_lowercase(Basis%NDindexq%sys_type)
      Call string_uppercase_TO_lowercase(Basis%NDindexb%sys_type)

      SELECT CASE (Basis%Basis_name)
      CASE ('dp')
        Allocate(Basis%tab_basis(nb_basis))
        Allocate(NDend_q(nb_basis))
        Allocate(NDend_b(nb_basis))
        Allocate(Basis%NDindexq%NDinit(nb_basis))
        Allocate(Basis%NDindexb%NDinit(nb_basis))

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

        Basis%nb = product(Basis%tab_basis(:)%nb)
        Basis%nq = product(Basis%tab_basis(:)%nq)

        DO i=1,NDim
          NDend_q(i)=Basis%tab_basis(i)%nq
          NDend_b(i)=Basis%tab_basis(i)%nb
          !NDend_l(i)=Basis%tab_basis(i)%nb
        END DO

        Call Init_NDindex(Basis%NDindexq,NDend_q,NDim)
        Call Init_NDindex(Basis%NDindexb,NDend_b,NDim)

      CASE ('smolyak')

        Basis%Smolyak%Ndim = nb_basis
        Allocate(NDend_l(Basis%Smolyak%Ndim))
        Allocate(Tab_ind(Basis%Smolyak%Ndim))
        Allocate(Tab_ind1(Basis%Smolyak%Ndim))
        !Allocate(Tab_li(Basis%Smolyak%Ndim))
        Allocate(Tags1(Basis%Smolyak%Ndim+1))

        NDend_l(:)=Basis%NDindexl%L
        Call Init_NDindex(Basis%NDindexl,NDend_l,Basis%Smolyak%Ndim)
        Call Init_tab_ind(Tab_ind,Basis%NDindexl)
        I=0
        Do
          I=I+1
          CALL increase_NDindex(Tab_ind,Basis%NDindexl,Endloop)
          Basis%Smolyak%Maxtermsmol = I
          write(out_unitp,*)I,Tab_ind
          IF (Endloop) exit
        END DO
        Deallocate(Tab_ind)
        write(out_unitp,*)"Basis%Smolyak%Maxtermsmol",Basis%Smolyak%Maxtermsmol

        Allocate(Basis%Smolyak%Li_smo(Basis%Smolyak%Maxtermsmol,Basis%Smolyak%Ndim))
        Allocate(Basis%Smolyak%DSmol(Basis%Smolyak%Maxtermsmol))

        Allocate(Basis%NDindexq_smol(Basis%Smolyak%Maxtermsmol))
        Allocate(Basis%NDindexb_smol(Basis%Smolyak%Maxtermsmol))

        Call Init_tab_ind(Tab_ind1,Basis%NDindexl)
        I=0
        Do
          I=I+1
          CALL increase_NDindex(Tab_ind1,Basis%NDindexl,Endloop)
          CALL Weight_D_smol (Basis%Smolyak%DSmol(I),Basis%NDindexl%L,sum(Tab_ind1),Basis%Smolyak%Ndim)

          Do inb=1,Basis%Smolyak%Ndim
            Basis%Smolyak%Li_smo(I,inb ) = tab_ind1(inb)
          END DO
          IF (Endloop) exit
        END DO
        Deallocate(Tab_ind1)

        Allocate(Basis%tab_Smolyak(Basis%Smolyak%Maxtermsmol))
        DO I=1,Basis%Smolyak%Maxtermsmol
          Allocate(Basis%tab_Smolyak(I)%tab_basis(Basis%Smolyak%Ndim))
        END DO

        Allocate(Basis%Smolyak%Nqlp(Basis%Smolyak%Maxtermsmol))
        Allocate(Basis%Smolyak%Nblp(Basis%Smolyak%Maxtermsmol))
        Allocate(Basis%Smolyak%Sum_Nqlp(Basis%Smolyak%Maxtermsmol))
        Allocate(Basis%Smolyak%Sum_Nblp(Basis%Smolyak%Maxtermsmol))
        Allocate(Basis%Smolyak1(Basis%Smolyak%Ndim))
        Allocate(Basis%Smolyak2(Basis%Smolyak%Maxtermsmol))
        Allocate(Basis%A_smol(Basis%Smolyak%Ndim))
        Allocate(Basis%B_smol(Basis%Smolyak%Ndim))


        DO J=1, Basis%Smolyak%Ndim
         Call Read_Basis_simple(name,nb_basis,LB,A_smol,B_smol,nb,nq,A,B,scaleQ,Q0,nio)

         Basis%A_smol(J)= A_smol
         Basis%B_smol(J)= B_smol

         Call Calc_n_smol_1(Basis,A_smol,B_smol,J)
         DO I=1,Basis%Smolyak%Maxtermsmol
           Basis%tab_Smolyak(I)%tab_basis(J)%nb_basis= nb_basis
           Basis%tab_Smolyak(I)%tab_basis(J)%A       = A
           Basis%tab_Smolyak(I)%tab_basis(J)%B       = B
           Basis%tab_Smolyak(I)%tab_basis(J)%scaleQ  = scaleQ
           Basis%tab_Smolyak(I)%tab_basis(J)%Q0      = Q0

           Basis%tab_Smolyak(I)%tab_basis(J)%nq = Calc_n_smol(A_smol,B_smol,Basis%Smolyak%Li_smo(I,J))
           Basis%tab_Smolyak(I)%tab_basis(J)%nb = Calc_n_smol(A_smol,B_smol,Basis%Smolyak%Li_smo(I,J))

           Basis%tab_Smolyak(I)%tab_basis(J)%Basis_name = trim(adjustl(name))
           Call Construct_basis_new (Basis%tab_Smolyak(I)%tab_basis(J))

         END DO
        END DO


        Allocate(Basis%Smolyak%NDend_smol_q(Basis%Smolyak%Ndim,Basis%Smolyak%Maxtermsmol))
        Allocate(Basis%Smolyak%NDend_smol_b(Basis%Smolyak%Ndim,Basis%Smolyak%Maxtermsmol))

        DO I=1, Basis%Smolyak%Maxtermsmol
        DO J=1, Basis%Smolyak%Ndim
          Basis%Smolyak%NDend_smol_q(J,I) = Basis%tab_Smolyak(I)%tab_basis(J)%nq
          Basis%Smolyak%NDend_smol_b(J,I) = Basis%tab_Smolyak(I)%tab_basis(J)%nb

        END DO
        END DO


        !Call Mapping_1 (Basis )


        DO I=1, Basis%Smolyak%Maxtermsmol
          Basis%NDindexq_smol(I)%sys_type = 'dp'
          Basis%NDindexb_smol(I)%sys_type = 'dp'
          Call Init_NDindex(Basis%NDindexq_smol(I),Basis%Smolyak%NDend_smol_q(:,I),Basis%Smolyak%Ndim)
          Call Init_NDindex(Basis%NDindexb_smol(I),Basis%Smolyak%NDend_smol_b(:,I),Basis%Smolyak%Ndim)
          Basis%Smolyak%Nblp(I)     = Product(Basis%tab_Smolyak(I)%tab_basis(:)%nb)
          Basis%Smolyak%Nqlp(I)     = Product(Basis%tab_Smolyak(I)%tab_basis(:)%nq)
          Basis%tab_Smolyak(I)%nb   = Product(Basis%tab_Smolyak(I)%tab_basis(:)%nb)
          Basis%tab_Smolyak(I)%nq   = Product(Basis%tab_Smolyak(I)%tab_basis(:)%nq)
        END DO

        DO I=1, Basis%Smolyak%Maxtermsmol
         IF(I==1) THEN
           Basis%Smolyak%Sum_Nqlp(I) =  Basis%Smolyak%Nqlp(I)
           Basis%Smolyak%Sum_Nblp(I) =  Basis%Smolyak%Nblp(I)
         ELSE
           Basis%Smolyak%Sum_Nqlp(I) = Basis%Smolyak%Sum_Nqlp(I-1) + Basis%Smolyak%Nqlp(I)
           Basis%Smolyak%Sum_Nblp(I) = Basis%Smolyak%Sum_Nblp(I-1) + Basis%Smolyak%Nblp(I)
         END IF
        END DO


        Basis%Smolyak%NSB       =  Sum(Basis%Smolyak%Nblp)
        Basis%Smolyak%NSQ       =  Sum(Basis%Smolyak%Nqlp)

        Call  Mapping_2G (Basis )
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
 END SUBROUTINE read_construct_basis

 FUNCTION Factoriel(n)
   Integer             :: Factoriel
   Integer, intent(in) :: n
   Integer             :: i

   IF(n==0) THEN
    Factoriel = 1
   ELSE IF(n<0) THEN
    Factoriel=1
   ELSE
    Factoriel=1
    Do i=1,n
     Factoriel = Factoriel*i
    END DO
   END IF
 END FUNCTION Factoriel

 FUNCTION Binomial(K,n)
   Integer                 :: Binomial
   Integer, intent(in)     :: K
   Integer, intent(in)     :: n
    Binomial = Factoriel(n)/(Factoriel(K)*Factoriel(n-K))
  !  Write(*,*)'n=',n,'k=',K,'Binomial=',Binomial
 END FUNCTION Binomial

 SUBROUTINE Weight_D_smol ( D,L,Som_l,NDIM)
  Integer,   intent(in)       :: NDIM
  Integer,   intent(in)       :: Som_l
  Integer,   intent(in)       :: L
  Integer,   intent(out)      :: D

  D = ((-1)**(L-Som_l))*Binomial(L-Som_l,NDIM-1)

END SUBROUTINE Weight_D_smol


 RECURSIVE SUBROUTINE Read_Basis_old (Basis,nio)
 USE UtilLib_m
   Logical,             parameter      :: debug = .true.
  !Logical,             parameter       ::debug = .false.
   TYPE(Basis_t),       intent(inout)  :: Basis
   Integer,             intent(in)     :: nio
   Integer, allocatable                :: NDend_q(:)
   Integer, allocatable                :: NDend_b(:)
   Integer, allocatable                :: NDend_l(:)
   Integer                             :: err_io,nb,nq,i,j,nb_basis,LB
   Character (len=Name_len)            :: name,sys_type
   Real(kind=Rk)                       :: A,B,scaleQ,Q0,d0,d2,X1,W1

   NAMELIST /basis_nD/ name,nb_basis,LB,nb,nq,A,B,scaleQ,Q0
   nb_basis  = 0
   nb        = 0
   nq        = 0
   LB        = 0
   A         = ZERO
   B         = ZERO
   Q0        = ZERO
   scaleQ    = ONE
   name      = '0'
   !sys_type  = 'dp'

   read(nio,nml=basis_nD,IOSTAT=err_io)
   write(out_unitp,nml=basis_nD)
   IF (err_io < 0) THEN
     write(out_unitp,basis_nD)
     write(out_unitp,*) ' ERROR in Read_Basis'
     write(out_unitp,*) ' while reading the namelist "basis_nD"'
     write(out_unitp,*) ' end of file or end of record'
     write(out_unitp,*) ' Probably, you forget a basis set ...'
     write(out_unitp,*) ' Check your data !!'
     STOP ' ERROR in Read_Basis: problems with the namelist.'
   END IF
   IF (err_io > 0) THEN
     write(out_unitp,basis_nD)
     write(out_unitp,*) ' ERROR in Read_Basis'
     write(out_unitp,*) ' while reading the namelist "basis_nD"'
     write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
     write(out_unitp,*) ' Check your data !!'
     STOP ' ERROR in Read_Basis: problems with the namelist.'
   END IF

   IF (nb_basis > 1) THEN

     Basis%Basis_name        = name
     Basis%NDindexq%sys_type = name
     Basis%NDindexb%sys_type = name
     Basis%NDindexl%sys_type = name
     Basis%NDindexb%L        = LB
     Basis%NDindexl%L        = LB
     Call string_uppercase_TO_lowercase(Basis%Basis_name)
     Call string_uppercase_TO_lowercase(Basis%NDindexl%sys_type)
     Call string_uppercase_TO_lowercase(Basis%NDindexq%sys_type)
     Call string_uppercase_TO_lowercase(Basis%NDindexb%sys_type)

     Allocate(Basis%tab_basis(nb_basis))
     Allocate(NDend_q(nb_basis))
     Allocate(NDend_b(nb_basis))
     Allocate(NDend_l(nb_basis))
     Allocate(Basis%NDindexq%NDinit(nb_basis))
     Allocate(Basis%NDindexb%NDinit(nb_basis))

     DO i=1,nb_basis

       Call Read_Basis_old(Basis%tab_basis(i),nio)
     END DO
     Basis%nb = product(Basis%tab_basis(:)%nb)
     Basis%nq = product(Basis%tab_basis(:)%nq)

     DO i=1,nb_basis
       NDend_q(i) = Basis%tab_basis(i)%nq
       NDend_b(i) = Basis%tab_basis(i)%nb
       NDend_l(i) = Basis%tab_basis(i)%nb
     END DO

     Call Init_NDindex(Basis%NDindexq,NDend_q,Size(Basis%tab_basis))
     Call Init_NDindex(Basis%NDindexb,NDend_b,Size(Basis%tab_basis))
     Call Init_NDindex(Basis%NDindexl,NDend_l,Size(Basis%tab_basis))
   ELSE
     Basis%nb_basis       = nb_basis
     Basis%nb             = nb
     Basis%nq             = nq
     Basis%A              = A
     Basis%B              = B
     Basis%Q0             = Q0
     Basis%scaleQ         = scaleQ
     Basis%Basis_name     = trim(adjustl(name))
     Call Construct_basis_new (Basis)
  END IF
END SUBROUTINE Read_Basis_old





 SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
 USE UtilLib_m
   TYPE(Basis_t),       intent(inout)  :: Basis
   Real(kind=Rk)                       :: dx
   Integer                             :: ib,iq,nb,nq

   nb = Basis%nb
   nq = Basis%nq
   dx = pi/nq
    ! grid and weight
   Basis%x = [(dx*(iq-HALF),iq=1,nq)]
   Basis%w = [(dx,iq=1,nq)]

   Allocate(Basis%d0gb(nq,nb))
   Allocate(Basis%d1gb(nq,nb,1))
   Allocate(Basis%d2gb(nq,nb,1,1))

   DO ib=1,nb
     Basis%d0gb(:,ib)     =          sin(Basis%x(:)*ib) / sqrt(pi*HALF)
     Basis%d1gb(:,ib,1)   =  ib    * cos(Basis%x(:)*ib) / sqrt(pi*HALF)
     Basis%d2gb(:,ib,1,1) = -ib**2 * Basis%d0gb(:,ib)
   END DO

   IF (nb == nq) THEN
     Basis%d0gb(:,nb)      = Basis%d0gb(:,nb)      / sqrt(TWO)
     Basis%d1gb(:,nb,:)    = Basis%d1gb(:,nb,:)    / sqrt(TWO)
     Basis%d2gb(:,nb,:,:)  = Basis%d2gb(:,nb,:,:)  / sqrt(TWO)
   END IF
 END SUBROUTINE Construct_Basis_Sin

 SUBROUTINE Construct_Basis_Ho(Basis) ! HO :
 USE UtilLib_m
  TYPE(Basis_t),       intent(inout)  :: Basis
  Integer                :: iq,ib

  Allocate(Basis%x(Basis%nq))
  Allocate(Basis%w(Basis%nq))

  Call hercom(Basis%nq, Basis%x(:), Basis%w(:))

  Allocate(Basis%d0gb(Basis%nq,Basis%nb))
  Allocate(Basis%d1gb(Basis%nq,Basis%nb,1))
  Allocate(Basis%d2gb(Basis%nq,Basis%nb,1,1))

  DO iq = 1, Basis%nq
    DO ib = 1, Basis%nb
      Call Construct_Basis_poly_Hermite_exp(Basis%x(iq),Basis%d0gb(iq,ib),&
      Basis%d1gb(iq,ib,1),Basis%d2gb(iq,ib,1,1), ib-1,.TRUE.)

    END DO
    !write(*,*)'Basis%d0gb(iq,ib)' ,Basis%d0gb(iq,:)
  END DO
 END SUBROUTINE Construct_Basis_Ho

 FUNCTION poly_Hermite(x,l)
  Implicit none
   Real(kind = Rk):: poly_Hermite
   Real(kind = Rk):: pl0,pl1,pl2,norme,x
   Integer        :: i,l

   IF ( l .LT. 0 ) THEN
    Write(out_unitp,*) 'Bad arguments in poly_hermite :'
    Write(out_unitp,*) ' l < 0 : ',l
    STOP
   END IF
    norme  =  sqrt(PI)
   IF (l .EQ. 0) THEN
    poly_Hermite = ONE/sqrt(norme)
   ELSE IF (l .EQ. 1) THEN
    norme = norme * TWO
    poly_Hermite = TWO * x/sqrt(norme)
   ELSE
    pl2 = ONE
    pl1 = TWO * x
    norme = norme * TWO
    DO i=2,l
     norme = norme * TWO * i
     pl0 = TWO*( x*pl1 - (i-1)*pl2 )
     pl2 = pl1
     pl1 = pl0
    END DO
    poly_Hermite = pl0/sqrt(norme)
   END IF
 END FUNCTION poly_Hermite

 FUNCTION gamma_perso(n)
  Implicit none
   Real(kind = Rk)  :: gamma_perso
   Real(kind = Rk)  :: a
   Integer          :: i,n
   IF (n .LE. 0) THEN
    write(out_unitp,*) 'ERROR: gamma( n<=0)',n
    STOP
   END IF
   a = ONE
   DO i = 1,n-1
    a = a * dble (i)
   END DO
   gamma_perso = a
 END FUNCTION gamma_perso

 SUBROUTINE herrec ( p2, dp2, p1, x, nq )
  Implicit none
  Integer        ::i
  Integer        :: nq
  Real(kind = Rk):: dp0,dp1,dp2,p0,p1,p2,x

  p1  = ONE
  dp1 = ZERO
  p2  = x
  dp2 = ONE

  DO i = 2, nq
    p0  = p1
    dp0 = dp1
    p1  = p2
    dp1 = dp2
    p2  = x * p1 - HALF * ( dble ( i ) - ONE ) * p0
    dp2 = x * dp1 + p1 - HALF * ( dble ( i ) - ONE ) * dp0
  END DO
 END SUBROUTINE herrec

  SUBROUTINE herroot ( x, nq, dp2, p1 )
  Implicit none
   Integer          :: i
   Integer          :: nq
   Real(kind = Rk),parameter :: eps = TEN**(-TWELVE) ! 1.0d-12
   Real(kind = Rk)  :: d,dp2,p1,p2,x

   DO i = 1, 10
     Call herrec ( p2, dp2, p1, x, nq )
     d = p2 / dp2
     x = x - d
     IF (ABS ( d ) .LE. eps * ( ABS ( x ) + ONE ) ) THEN
      RETURN
     END IF
   END DO
 END SUBROUTINE herroot

 SUBROUTINE hercom (nq,xp,w)
 Implicit none
   Integer        :: i,nq
   Real(kind = Rk):: cc,dp2,p1,s,temp,x
   Real(kind = Rk):: w(nq),xp(nq)

   CC = 1.7724538509_Rk * gamma_perso(nq ) / ( TWO**( nq-1) )

   S = ( TWO * dble (Real(nq,Kind=Rk) ) + ONE )**( SIXTH )

   DO i = 1, ( nq + 1 ) / 2
     IF ( i .EQ. 1 ) THEN
      x = s**3 - 1.85575_Rk / s
     ELSE IF ( i .EQ. 2 ) THEN
      x = x - 1.14_Rk * ( ( dble ( nq ) )**0.426_Rk ) / x
     ELSE IF ( i .EQ. 3 ) THEN
      x = 1.86_Rk * x - 0.86_Rk * xp(1)
     ELSE IF ( i .EQ. 4 ) THEN
      x = 1.91_Rk * x - 0.91_Rk * xp(2)
     ELSE
      x = TWO * x - xp(i-2)
     END IF
     Call herroot ( x,  nq, dp2, p1 )
     xp(i) = x
     W(i) = cc / dp2 / p1
     xp( nq-i+1) = - x
     w( nq-i+1) = w(i)
   END DO
   DO i = 1,  nq/2
     temp = xp(i)
     xp(i) = xp( nq+1-i)
     xp( nq+1-i) = temp
   END DO
   DO i = 1, nq
     w(i) = w(i)*exp(xp(i)*xp(i))
   END DO
 END SUBROUTINE hercom

 SUBROUTINE Construct_Basis_poly_Hermite_exp(x,d0gb,d1gb,d2gb,l,deriv)
  Logical        :: deriv
  Integer        :: l
  Real(kind = RK):: pexp,x,d0gb,d1gb,d2gb

  IF (deriv) THEN
    d0gb = poly_Hermite( x,l)
    IF (l .EQ. 0) THEN
     d1gb     = ZERO
     d2gb     = ZERO
    ELSE IF (l .EQ. 1) THEN
     d1gb = sqrt(TWO)*poly_Hermite( x,0)
     d2gb = ZERO
    ELSE IF (l .EQ. 2) THEN
     d1gb = sqrt(TWO*l) * poly_Hermite( x,l-1)
     d2gb = TWO*( x*d1gb-d0gb *l)
    ELSE
     d1gb = sqrt(TWO*l) * poly_Hermite( x,l-1)
     d2gb = TWO*( x* d1gb-d0gb*l)
    END IF
    pexp = exp(- HALF* x* x)
    d2gb = (d2gb-TWO*x*d1gb+( x* x-ONE)*d0gb)*pexp
    d1gb = (d1gb- x*d0gb)*pexp
    d0gb = d0gb*pexp
  ELSE
    d0gb = poly_Hermite(x ,l)*exp(-HALF* x* x)
    d1gb = ZERO
    d2gb = ZERO
  END IF
 END SUBROUTINE Construct_Basis_poly_Hermite_exp

 SUBROUTINE CheckOrtho_Basis(Basis,nderiv)
 USE UtilLib_m
   Logical,                 parameter   :: debug = .true.
   !Logical,                parameter   ::debug = .false.
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
    IF (nderiv > -1) Call Write_RMat(S,out_unitp,5,name_info='S')
    Sii = ZERO
    Sij = ZERO
    DO ib=1,Basis%nb
      IF (abs(S(ib,ib)-ONE) > Sii) Sii = abs(S(ib,ib)-ONE)
      S(ib,ib) = ZERO
    END DO

    Sij = maxval(S)
    write(out_unitp,*) 'Sii,Sij',Sii,Sij

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
    Logical,          parameter     :: debug = .false.
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
   real(kind=Rk), ALLOCATABLE   :: d0bgw(:,:)
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
    logical          , parameter             :: debug = .true.
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
    Integer                         :: ib
    !Logical,          parameter    :: debug = .true.
    Logical,         parameter    ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Calc_dngg_grid'
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
      write(out_unitp,*) 'END Calc_dngg_grid'
      flush(out_unitp)
    END IF
    deAllocate(d0bgw)
  END SUBROUTINE Calc_dngg_grid

  SUBROUTINE Test_Passage(Basis)
  USE UtilLib_m
    TYPE(Basis_t),    intent(in)    :: Basis
    Logical,          parameter    :: debug = .true.
    !Logical,         parameter    ::debug = .false.
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

    x0 = Basis%Q0
    sx = Basis%scaleQ
    IF (abs(sx) > ONETENTH**6 .AND. Basis_IS_Allocated(Basis)) THEN

      Basis%x(:) = x0 + Basis%x(:) / sx
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
