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
  PUBLIC :: Basis_t,Read_Basis,Basis_IS_Allocated,BasisTOGrid_Basis,GridTOBasis_Basis
  PUBLIC :: Test_Passage,Calc_dngg_grid,Basis_IS_Allocatedtot,write_basis
  PUBLIC :: BasisTOGrid_Basis_rapide,GridTOBasis_Basis_rapide
  PUBLIC :: BasisTOGrid_Basis_rapide1,GridTOBasis_Basis_rapide1

  TYPE :: Basis_t
    ! nb_basis is number of bases used
    Integer                      :: nb_basis   = 0
    ! nb  is number of basic elements according to each degree of freedom
    Integer                      :: nb         = 0
    ! nq  is number of grid points along each degree of freedom
    Integer                      :: nq         = 0
    ! Basis_name is the name of the base  used
    Character(len=:),allocatable :: Basis_name
    ! x(:) are the positions of the grid points along each degree of freedom
    Real(kind=Rk),   allocatable :: x(:)
    ! w(:) are the weights of the grid points along each degree of freedom
    Real(kind=Rk),   allocatable :: w(:)
    Real(kind=Rk),   allocatable :: d0gb(:,:)      ! basis functions d0gb(nq,nb)
    Real(kind=Rk),   allocatable :: d1gb(:,:,:)    ! basis functions d2gb(nq,nb,1)
    Real(kind=Rk),   allocatable :: d1gg(:,:,:)    ! basis functions d2gg(nq,nq,1)
    Real(kind=Rk),   allocatable :: d2gb(:,:,:,:)  ! basis functions d2gb(nq,nb,1,1)
    Real(kind=Rk),   allocatable :: d2gg(:,:,:,:)  ! basis functions d2gg(nq,nq,1,1)
    TYPE(NDindex_t)              :: NDindexq
    TYPE(NDindex_t)              :: NDindexb
    TYPE (Basis_t),  allocatable :: tab_basis(:)
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

  RECURSIVE SUBROUTINE Read_Basis(Basis,nio)
  USE UtilLib_m
    Logical,             parameter     :: debug = .true.
   !Logical,             parameter     ::debug = .false.
    TYPE(Basis_t),       intent(inout)  :: Basis
    Integer,             intent(in)     :: nio
    Integer, allocatable                :: NDend_q(:)
    Integer, allocatable                :: NDend_b(:)
    Integer                             :: err_io,nb,nq,i,j,nb_basis
    Character (len=Name_len)            :: name
    Real(kind=Rk)                       :: A,B,scaleQ,Q0,d0,d2,X1,W1

    NAMELIST /basis_nD/ name,nb_basis,nb,nq,A,B,scaleQ,Q0
    nb_basis  = 0
    nb        = 0
    nq        = 0
    A         = ZERO
    B         = ZERO
    Q0        = ZERO
    scaleQ    = ONE
    name      = '0'

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
      Basis%Basis_name     = 'dp'

      Call string_uppercase_TO_lowercase(Basis%Basis_name)

      Allocate(Basis%tab_basis(nb_basis))
      Allocate(NDend_q(nb_basis))
      Allocate(NDend_b(nb_basis))

      DO i=1,nb_basis
        Call Read_Basis(Basis%tab_basis(i),nio)
      END DO
      Basis%nb = product(Basis%tab_basis(:)%nb)
      Basis%nq = product(Basis%tab_basis(:)%nq)

      DO i=1,nb_basis
        NDend_q(i)=Basis%tab_basis(i)%nq
        NDend_b(i)=Basis%tab_basis(i)%nb
      END DO

      Call Init_NDindex(Basis%NDindexq,NDend_q,nb_basis)
      Call Init_NDindex(Basis%NDindexb,NDend_b,nb_basis)

    ELSE
      Basis%nb_basis  = nb_basis
      Basis%nb        = nb
      Basis%nq        = nq
      Basis%Basis_name     = trim(adjustl(name))
      Call string_uppercase_TO_lowercase(Basis%Basis_name)

      SELECT CASE (Basis%Basis_name)
      CASE ('boxab')
      Call Construct_Basis_Sin(Basis)
      Q0      = A
      scaleQ  = pi/(B-A)
      CASE ('herm','ho')
      Call Construct_Basis_Ho(Basis)
      CASE default
      STOP 'ERROR in Read_Basis: no default basis.'
    END SELECT
      Call Scale_Basis(Basis,Q0,scaleQ)
      Call Calc_dngg_grid(Basis)
      Call CheckOrtho_Basis(Basis,nderiv=2)
   END IF
 END SUBROUTINE Read_Basis

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

!!!!modife!!12:01:2023!
  SUBROUTINE Vect_change(BG1,BG2)
    Real (kind=Rk),  intent(inout)       :: BG1(:)
    Real (kind=Rk),  intent(inout)       :: BG2(:)
        BG1(:) = BG2(:)
        BG2(:) = ZERO

  END SUBROUTINE Vect_change

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
    integer                                  :: i1,i3,iq,ib

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

  SUBROUTINE Calc_indice( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
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

      Iq1(1) = 1
      Ib1(1) = 1

      Iq2(1) =  Basis%tab_basis(1)%nq
      Ib2(1) =  Basis%tab_basis(1)%nb

      Iq3(1) =  Product(Basis%tab_basis(2:Ndim)%nq)
      Ib3(1) =  Product(Basis%tab_basis(2:Ndim)%nb)


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

  END SUBROUTINE Calc_indice




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
    Real (kind=Rk) ,allocatable,target      :: GB1(:),GB2(:)
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
      GBB(1:1,1:Basis%tab_basis(1)%nq,1:Basis%tab_basis(2)%nb*Basis%tab_basis(3)%nb) =>  GBB1

      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(1))
!

      Allocate(GGB2(Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nq*Basis%tab_basis(3)%nb))

      BBB(1:Basis%tab_basis(1)%nq,1:Basis%tab_basis(2)%nb,1:Basis%tab_basis(3)%nb) =>  GBB1
      GBB(1:Basis%tab_basis(1)%nq,1:Basis%tab_basis(2)%nq,1:Basis%tab_basis(3)%nb) =>  GGB2

      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(2))

      Deallocate(GBB1)


      BBB(1:Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nq,1:Basis%tab_basis(3)%nb,1:1) => GGB2
      GBB(1:Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nq,1:Basis%tab_basis(3)%nq,1:1 )   => G

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
      Call Calc_indice( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)

      Allocate(GBB1(iq1(1)*iq2(1)*ib3(1)))

      BBB(1:iq1(1),1:ib2(1),1:ib3(1)) => B
      GBB(1:iq1(1),1:iq2(1),1:ib3(1)) => GBB1
      !write(out_unitp,*) 'OK0'
      !write(out_unitp,*)ib1(1), ib2(1),ib3(1)

      !write(out_unitp,*)iq1(1), iq2(1),iq3(1)
      Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(1))

      DO inb = 2,Ndim-1

        Allocate(GGB2(iq1(inb)*iq2(inb)*ib3(inb)))

        BBB(1:iq1(Inb),1:ib2(inb),1:ib3(inb)) => GBB1

        GBB(1:iq1(inb),1:iq2(inb),1:ib3(inb)) => GGB2

        Call BasisTOGrid_1D(GBB,BBB,Basis%tab_basis(inb))

         GBB1=GGB2
        !Deallocate(GBB1)

        !Allocate(GBB1(iq1(inb)*iq2(inb)*ib3(inb)))

        !Call Vect_change(GBB1,GGB2)

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
!!!!fin modife!!12:01:2023!

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

  !!!!modife!!12:01:2023!

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
      Call Calc_indice( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Basis)
      Allocate(BGG1(Iq1(1)*Iq2(1)*Iq3(1)))
      BGG1(:) = ZERO
      GGG( 1:Iq1(1),1:Iq2(1),1:Iq3(1))   => BGG1
      GGB(1:Iq1(1),1:Iq2(1),1:Ib3(1))    => G

      Call GridTOBasis_1D(GGG,GGB,Basis%tab_basis(1))

      DO inb = 2,Ndim-1
        Allocate(BGG2(Ib1(inb)*Ib2(inb)*Iq3(inb)))
        BGG2(:) = ZERO
        GGG( 1:Iq1(inb),1:Iq2(inb),1:Ib3(inb))    => BGG2
        GGB( 1:Iq1(inb),1:Iq2(inb),1:Ib3(inb))    => BGG1
        Call GridTOBasis_1D(GGG,GGB,Basis%tab_basis(inb))
        Deallocate(BGG1)
        Allocate(BGG1(Ib1(inb)*Ib2(inb)*Iq3(inb)))
        Call Vect_change(BGG2,BGG1)
        Deallocate(BGG2)
      END DO

      B(:) = ZERO

      GGB(1:Iq1(Ndim),1:Ib2(Ndim),1:Ib3(Ndim)) => BGG1
      GGG(1:Ib1(Ndim),1:Ib2(Ndim),1:Ib3(Ndim)) => B
      Call GridTOBasis_1D(GGG,GGB,Basis%tab_basis(Ndim))

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
!!!!fin modife!!12:01:2023!
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

SUBROUTINE Scale_Basis(Basis,x0,sx)
USE UtilLib_m
    TYPE(Basis_t),       intent(inout)  :: Basis
    Real(kind=Rk),       intent(in)     :: x0,sx

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
