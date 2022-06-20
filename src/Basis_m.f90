MODULE Basis_m
  USE NumParameters_m
  USE NDindex_m
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Basis_t,Read_Basis,Basis_IS_allocated,BasisTOGrid_Basis,GridTOBasis_Basis
  PUBLIC :: Test_Passage,Calc_dngg_grid,Basis_IS_allocatedtot,write_basis

  TYPE :: Basis_t
    integer                      :: nb_basis   = 0
    integer                      :: nb         = 0
    integer                      :: nq         = 0
    character(len=:),allocatable :: Basis_name
    real(kind=Rk),   allocatable :: x(:)
    real(kind=Rk),   allocatable :: w(:)
    real(kind=Rk),   allocatable :: d0gb(:,:)      ! basis functions d0gb(nq,nb)
    real(kind=Rk),   allocatable :: d1gb(:,:,:)    ! basis functions d2gb(nq,nb,1)
    real(kind=Rk),   allocatable :: d1gg(:,:,:)    ! basis functions d2gg(nq,nq,1)
    real(kind=Rk),   allocatable :: d2gb(:,:,:,:)  ! basis functions d2gb(nq,nb,1,1)
    real(kind=Rk),   allocatable :: d2gg(:,:,:,:)  ! basis functions d2gg(nq,nq,1,1)
    TYPE(NDindex_t)              :: NDindexq
    TYPE(NDindex_t)              :: NDindexb
    TYPE (Basis_t),  allocatable :: tab_basis(:)
  END TYPE Basis_t

CONTAINS
RECURSIVE FUNCTION Basis_IS_allocated(Basis) RESULT(alloc)

    TYPE(Basis_t),   intent(in)  :: Basis
    logical                      :: alloc
    integer                      :: i

    alloc = allocated(Basis%tab_basis)
    IF ( allocated(Basis%tab_basis)) THEN
      Do i=1,size(Basis%tab_basis)
        alloc  = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
      END DO
    ELSE
      alloc =             allocated(Basis%x)
      alloc = alloc .AND. allocated(Basis%w)
      alloc = alloc .AND. allocated(Basis%d0gb)
      alloc = alloc .AND. allocated(Basis%d1gb)
      alloc = alloc .AND. allocated(Basis%d2gb)
    END IF
  END FUNCTION Basis_IS_allocated

  RECURSIVE FUNCTION Basis_IS_allocatedtot(Basis) RESULT(alloc)

      TYPE(Basis_t),   intent(in)  :: Basis
      logical                      :: alloc
      integer                      :: i

      alloc = allocated(Basis%tab_basis)
      IF ( allocated(Basis%tab_basis)) THEN
        Do i=1,size(Basis%tab_basis)
          alloc  = alloc .and. Basis_IS_allocated(Basis%tab_basis(i))
        END DO
      ELSE
        alloc =             allocated(Basis%x)
        alloc = alloc .AND. allocated(Basis%w)
        alloc = alloc .AND. allocated(Basis%d0gb)
        alloc = alloc .AND. allocated(Basis%d1gb)
        alloc = alloc .AND. allocated(Basis%d2gb)
        alloc = alloc .AND. allocated(Basis%d1gg)
        alloc = alloc .AND. allocated(Basis%d2gg)
      END IF

    END FUNCTION Basis_IS_allocatedtot

  RECURSIVE SUBROUTINE Write_Basis(Basis)
  USE UtilLib_m

    TYPE(Basis_t),       intent(in)  :: Basis
    integer                          :: i

    write(out_unitp,*) '-------------------------------------------------'
    write(out_unitp,*) 'Write_Basis'
    write(out_unitp,*) 'nb,nq',Basis%nb,Basis%nq


      IF (.NOT.allocated(Basis%x)) THEN
       write(out_unitp,*)' Basis table x is not allocated.'
      ELSE
        CALL Write_RVec(Basis%x,out_unitp,5,name_info='x')
      END IF
      write(out_unitp,*)
      IF (.NOT.allocated(Basis%W)) THEN
        write(out_unitp,*)' Basis table w is not allocated.'
      ELSE
        CALL Write_RVec(Basis%w,out_unitp,5,name_info='w')
      END IF
      write(out_unitp,*)
      IF (.NOT.allocated(Basis%d0gb)) THEN
        write(out_unitp,*)' Basis table d0gb is not allocated.'
      ELSE
        CALL Write_RMat(Basis%d0gb,out_unitp,5,name_info='d0gb')
      END IF
      write(out_unitp,*)
      IF (.NOT.allocated(Basis%d1gb)) THEN
        write(out_unitp,*)' Basis table d1gb is not allocated.'
      ELSE
        CALL Write_RMat(Basis%d1gb(:,:,1),out_unitp,5,name_info='d1gb')
      END IF
      write(out_unitp,*)
      IF (.NOT.allocated(Basis%d1gg)) THEN
        write(out_unitp,*)' Basis table d1gb is not allocated.'
      ELSE
        CALL Write_RMat(Basis%d1gg(:,:,1),out_unitp,5,name_info='d1gg')
      END IF
      write(out_unitp,*)
      IF (.NOT.allocated(Basis%d2gb)) THEN
        write(out_unitp,*)' Basis table d1gb is not allocated.'
      ELSE
        CALL Write_RMat(Basis%d2gb(:,:,1,1),out_unitp,5,name_info='d2gb')
      END IF
      write(out_unitp,*)
      IF (.NOT.allocated(Basis%d2gg)) THEN
        write(out_unitp,*)' Basis table d2gg is not allocated.'
      ELSE
        CALL Write_RMat(Basis%d2gg(:,:,1,1),out_unitp,5,name_info='d2gg')
      END IF

    !  write(out_unitp,*) 'nb_basis',Basis%nb_basis
    IF (allocated(Basis%tab_basis)) THEN
      DO i=1,size(Basis%tab_basis)
        CALL Write_Basis(Basis%tab_basis(i))
      END DO
    END IF
    write(out_unitp,*) '-------------------------------------------------'

  END SUBROUTINE Write_Basis

  RECURSIVE SUBROUTINE Read_Basis(Basis,nio)
  USE UtilLib_m
    logical,             parameter      :: debug = .true.
   !logical,             parameter      ::debug = .false.
    TYPE(Basis_t),       intent(inout)  :: Basis
    integer,             intent(in)     :: nio
    integer, allocatable                :: NDend_q(:)
    integer, allocatable                :: NDend_b(:)
    integer                             :: err_io,nb,nq,i,j,nb_basis
    character (len=Name_len)            :: name
    real(kind=Rk)                       :: A,B,scaleQ,Q0,d0,d2,X1,W1

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
      Basis%Basis_name     = 'Dp'

      CALL string_uppercase_TO_lowercase(Basis%Basis_name)

      allocate(Basis%tab_basis(nb_basis))
      allocate(NDend_q(nb_basis))
      allocate(NDend_b(nb_basis))

      DO i=1,nb_basis
        CALL Read_Basis(Basis%tab_basis(i),nio)
      END DO
      Basis%nb = product(Basis%tab_basis(:)%nb)
      Basis%nq = product(Basis%tab_basis(:)%nq)

      DO i=1,nb_basis
        NDend_q(i)=Basis%tab_basis(i)%nq
        NDend_b(i)=Basis%tab_basis(i)%nb
      END DO

      CALL Init_NDindex(Basis%NDindexq,NDend_q,nb_basis)
      CALL Init_NDindex(Basis%NDindexb,NDend_b,nb_basis)

    ELSE
      Basis%nb_basis  = nb_basis
      Basis%nb        = nb
      Basis%nq        = nq
      Basis%Basis_name     = trim(adjustl(name))
      CALL string_uppercase_TO_lowercase(Basis%Basis_name)

    SELECT CASE (Basis%Basis_name)
      CASE ('boxab')
      CALL Construct_Basis_Sin(Basis)
      Q0      = A
      scaleQ  = pi/(B-A)
      CASE ('herm','ho')
      CALL Construct_Basis_Ho(Basis)
      CASE default
      STOP 'ERROR in Read_Basis: no default basis.'
    END SELECT
      CALL Scale_Basis(Basis,Q0,scaleQ)
      CALL Calc_dngg_grid(Basis)
      CALL CheckOrtho_Basis(Basis,nderiv=2)
   END IF
 END SUBROUTINE Read_Basis

 SUBROUTINE Construct_Basis_Sin(Basis) ! sin : boxAB with A=0 and B=pi
 USE UtilLib_m
   TYPE(Basis_t),       intent(inout)  :: Basis
   real(kind=Rk)                       :: dx
   integer                             :: ib,iq,nb,nq

   nb = Basis%nb
   nq = Basis%nq
   dx = pi/nq
    ! grid and weight
   Basis%x = [(dx*(iq-HALF),iq=1,nq)]
   Basis%w = [(dx,iq=1,nq)]

   allocate(Basis%d0gb(nq,nb))
   allocate(Basis%d1gb(nq,nb,1))
   allocate(Basis%d2gb(nq,nb,1,1))

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
  integer                :: iq,ib

  allocate(Basis%x(Basis%nq))
  allocate(Basis%w(Basis%nq))

  call hercom(Basis%nq, Basis%x(:), Basis%w(:))

  allocate(Basis%d0gb(Basis%nq,Basis%nb))
  allocate(Basis%d1gb(Basis%nq,Basis%nb,1))
  allocate(Basis%d2gb(Basis%nq,Basis%nb,1,1))

  DO iq = 1, Basis%nq
    DO ib = 1, Basis%nb
      CALL Construct_Basis_poly_Hermite_exp(Basis%x(iq),Basis%d0gb(iq,ib),&
      Basis%d1gb(iq,ib,1),Basis%d2gb(iq,ib,1,1), ib-1,.TRUE.)
    END DO
  END DO
 END SUBROUTINE Construct_Basis_Ho

 FUNCTION poly_Hermite(x,l)
  Implicit none
   real(kind = Rk):: poly_Hermite
   real(kind = Rk):: pl0,pl1,pl2,norme,x
   integer        :: i,l

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
   real(kind = Rk)  :: gamma_perso
   real(kind = Rk)  :: a
   integer          :: i,n
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
  integer        ::i
  integer        :: nq
  real(kind = Rk):: dp0,dp1,dp2,p0,p1,p2,x

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
   integer          :: i
   integer          :: nq
   real(kind = Rk),parameter  :: eps = TEN**(-TWELVE) ! 1.0d-12
   real(kind = Rk)  :: d,dp2,p1,p2,x

   DO i = 1, 10
     CALL herrec ( p2, dp2, p1, x, nq )
     d = p2 / dp2
     x = x - d
     IF (ABS ( d ) .LE. eps * ( ABS ( x ) + ONE ) ) THEN
      RETURN
     END IF
   END DO
 END SUBROUTINE herroot

 SUBROUTINE hercom (nq,xp,w)
 Implicit none
   integer        :: i,nq
   real(kind = Rk):: cc,dp2,p1,s,temp,x
   real(kind = Rk):: w(nq),xp(nq)

   CC = 1.7724538509_Rk * gamma_perso(nq ) / ( TWO**( nq-1) )

   S = ( TWO * dble (real(nq,Kind=Rk) ) + ONE )**( SIXTH )

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
     CALL herroot ( x,  nq, dp2, p1 )
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
  logical        :: deriv
  integer        :: l
  real(kind = RK):: pexp,x,d0gb,d1gb,d2gb

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
   logical,                 parameter    :: debug = .true.
   !logical,                parameter    ::debug = .false.
   TYPE(Basis_t),           intent(in)   :: Basis
   integer,                 intent(in)   :: nderiv
   integer                               :: ib
   real(kind=Rk), ALLOCATABLE            :: S(:,:)
   real(kind=Rk), ALLOCATABLE            :: d0bgw(:,:)
   real(kind=Rk)                         :: Sii,Sij

   IF (Basis_IS_allocated(Basis)) THEN
    d0bgw = transpose(Basis%d0gb)
    DO ib=1,Basis%nb
      d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
    END DO
    S = matmul(d0bgw,Basis%d0gb)
    IF (nderiv > -1) CALL Write_RMat(S,out_unitp,5,name_info='S')
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
      !CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>',Rformat='e13.4')
      CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
    END IF

    IF (nderiv > 1) THEN
      write(out_unitp,*)
      S = matmul(d0bgw,Basis%d2gb(:,:,1,1))
      !CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d2b>',Rformat='e13.4')
      CALL Write_RMat(S,out_unitp,5,name_info='<d0b|d1b>')
    END IF
  ELSE
    write(out_unitp,*) ' WARNNING in CheckOrtho_Basis'
    write(out_unitp,*) ' the basis is not allocated.'
  END IF
  END SUBROUTINE CheckOrtho_Basis

  SUBROUTINE BasisTOGrid_Basis(G,B,Basis)
  USE UtilLib_m
  USE NDindex_m
    !logical,           parameter     :: debug = .true.
    logical,          parameter      :: debug = .false.
    TYPE(Basis_t),     intent(in)    :: Basis
    real(kind=Rk),     intent(in)    :: B(:)
    real(kind=Rk),     intent(inout) :: G(:)
    real(kind=Rk)                    :: W
    logical                          :: Endloop_q
    logical                          :: Endloop_b
    integer,         allocatable     :: tab_iq(:)
    integer,         allocatable     :: tab_ib(:)
    integer                          :: ib,iq,nq,nb,inb
    integer                          :: jb,jb1,jb2

    IF(debug)THEN
     write(out_unitp,*) 'BEGINNING BasisTOGrid_Basis'
     write(out_unitp,*) 'intent(in) :: B(:)',B
     CALL Write_Basis(Basis)
     flush(out_unitp)
    END IF

    IF(.NOT. Basis_IS_allocated(Basis))THEN
     write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
     write(out_unitp,*) " the basis is not allocated."
     STOP "ERROR BasisTOGrid_Basis: the basis is not allocated."
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

    IF(allocated(Basis%tab_basis)) THEN
     Allocate(Tab_ib(size(Basis%tab_basis)))
     Allocate(Tab_iq(size(Basis%tab_basis)))
     Call Init_tab_ind(Tab_iq,Basis%NDindexq)
     Iq=0
     DO
      Iq=Iq+1
      CALL increase_NDindex(Tab_iq,Basis%NDindexq,Endloop_q)
      IF (Endloop_q) exit
      G(iq)=ZERO
      Call Init_tab_ind(Tab_ib,Basis%NDindexb)
      Ib=0
      DO
       Ib=Ib+1
       CALL increase_NDindex(Tab_ib,Basis%NDindexb,Endloop_b)
       IF (Endloop_b) exit
       W=ONE
       DO inb=1,size(Basis%tab_basis)
        W=W* Basis%tab_basis(inb)%d0gb(tab_iq(inb),tab_ib(inb))
       END DO
       G(iq) =G(iq)+W*B(ib)
      END DO
     END DO
     Deallocate(Tab_ib)
     Deallocate(Tab_iq)
    ELSE
     DO iq=1,Basis%nq
      G(iq)=ZERO
      DO ib=1,Basis%nb
       G(iq)= G(iq)+Basis%d0gb(iq,ib)*B(ib)
      END DO
     END DO
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'intent(OUTIN) :: G(:)',G
      write(out_unitp,*) 'END BasisTOGrid_Basis'
      flush(out_unitp)
    END IF
  END SUBROUTINE BasisTOGrid_Basis

  SUBROUTINE GridTOBasis_Basis(B,G,Basis)
  USE UtilLib_m
    !logical,          parameter     :: debug = .true.
    logical,         parameter      :: debug = .false.
    TYPE(Basis_t),    intent(in)    :: Basis
    real(kind=Rk),    intent(in)    :: G(:)
    real(kind=Rk),    intent(inout) :: B(:)
    logical                         :: Endloop_q
    logical                         :: Endloop_b
    real(kind=Rk)                   :: WT,W
    integer,        allocatable     :: tab_iq(:)
    integer,        allocatable     :: tab_ib(:)
    integer                         :: ib,iq,iq1,iq2,inb
    integer                         :: jb,ib1,ib2,jb1,jb2

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING GridTOBasis_Basis'
      write(out_unitp,*) 'intent(in) :: G(:)',G
      !CALL Write_Basis(Basis)
      flush(out_unitp)
    END IF

    IF (.NOT. Basis_IS_allocated(Basis)) THEN
      write(out_unitp,*) ' ERROR in BasisTOGrid_Basis'
      write(out_unitp,*) " the basis is not allocated."
      STOP "ERROR BasisTOGrid_Basis: the basis is not allocated."
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

    IF(allocated(Basis%tab_basis)) THEN
     Allocate(Tab_ib(size(Basis%tab_basis)))
     Allocate(Tab_iq(size(Basis%tab_basis)))
     Call Init_tab_ind(Tab_ib,Basis%NDindexb)
     Ib=0
     DO
      Ib=Ib+1
      CALL increase_NDindex(Tab_ib,Basis%NDindexb,Endloop_b)
      IF (Endloop_b) exit
      B(ib)=ZERO
      Call Init_tab_ind(Tab_iq,Basis%NDindexq)
      Iq=0
      DO
        Iq=Iq+1
        CALL increase_NDindex(Tab_iq,Basis%NDindexq,Endloop_q)
        IF (Endloop_q) exit
        WT=1
        W=1
        DO inb=1,size(Basis%tab_basis)
         WT= WT *Basis%tab_basis(inb)%w(tab_iq(inb))
         W=W*Basis%tab_basis(inb)%d0gb(tab_iq(inb),tab_ib(inb))
        END DO
        B(ib)=B(ib)+W*WT*G(iq)
      END DO
     END DO
     Deallocate(Tab_ib)
     Deallocate(Tab_iq)
    ELSE
     DO ib=1,Basis%nb
       B(ib)=ZERO
       DO iq=1,Basis%nq
         B(ib)=B(ib)+Basis%d0gb(iq,ib)*Basis%w(iq)*G(iq)
       END DO
     END DO
    END IF

    IF (debug) THEN
     ! write(out_unitp,*) 'intent(OUTIN) :: B(:)',B
     write(out_unitp,*) 'END GridTOBasis_Basis'
     flush(out_unitp)
    END IF
  END SUBROUTINE GridTOBasis_Basis

  SUBROUTINE Calc_dngg_grid(Basis)
  USE UtilLib_m
    TYPE(Basis_t), intent(inout)    :: Basis
    real(kind=Rk), allocatable      :: d0bgw(:,:)
    integer                         :: ib
    !logical,          parameter     :: debug = .true.
    logical,         parameter     ::debug = .false.

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Calc_dngg_grid'
      CALL Write_Basis(Basis)
      flush(out_unitp)
    END IF

    allocate(Basis%d1gg(Basis%nq,Basis%nq,1))
    allocate(Basis%d2gg(Basis%nq,Basis%nq,1,1))

    d0bgw = transpose(Basis%d0gb)
    DO ib=1,Basis%nb
       d0bgw(ib,:) = d0bgw(ib,:) * Basis%w(:)
    END DO

    IF (debug) THEN
      CALL Write_RMat(d0bgw(:,:),out_unitp,5,name_info='d0bgw')
      write(out_unitp,*)
    END IF

    Basis%d1gg(:,:,1)   = matmul(Basis%d1gb(:,:,1),d0bgw)
    Basis%d2gg(:,:,1,1) = matmul(Basis%d2gb(:,:,1,1),d0bgw)

    !CALL Write_RMat(Basis%d1gg(:,:,1),out_unitp,5,name_info='d1gg')
    !write(out_unitp,*)
    !CALL Write_RMat(Basis%d2gg(:,:,1,1),out_unitp,5,name_info='d2gg')

    IF (debug) THEN
      CALL Write_Basis(Basis)
      write(out_unitp,*) 'END Calc_dngg_grid'
      flush(out_unitp)
    END IF
    deallocate(d0bgw)
  END SUBROUTINE Calc_dngg_grid

  SUBROUTINE Test_Passage(Basis)
  USE UtilLib_m
    TYPE(Basis_t),    intent(in)    :: Basis
    logical,          parameter     :: debug = .true.
    !logical,         parameter     ::debug = .false.
    real(kind=Rk),    allocatable   :: G1(:),B1(:)
    real(kind=Rk),    allocatable   :: G2(:),B2(:)
    real(kind=Rk),    allocatable   :: B(:),Delta_g(:),Delta_b(:)
    real(kind=Rk),    parameter     :: eps = TEN**(-TEN)
    integer                         :: iq

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING Test_Passage'
      CALL Write_Basis(Basis)
      flush(out_unitp)
    END IF

    allocate(B(Basis%nb))
    allocate(Delta_g(Basis%nq))
    allocate(G1(Basis%nq))
    allocate(G2(Basis%nq))

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
      allocate(Delta_b(Basis%nb))
      allocate(B1(Basis%nb))
      allocate(B2(Basis%nb))
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
      CALL Write_Basis(Basis)
      flush(out_unitp)
    END IF
END SUBROUTINE Test_Passage

SUBROUTINE Scale_Basis(Basis,x0,sx)
USE UtilLib_m
    TYPE(Basis_t),       intent(inout)  :: Basis
    real(kind=Rk),       intent(in)     :: x0,sx

    IF (abs(sx) > ONETENTH**6 .AND. Basis_IS_allocated(Basis)) THEN

      Basis%x(:) = x0 + Basis%x(:) / sx
      Basis%w(:) =      Basis%w(:) / sx

      Basis%d0gb(:,:)     = Basis%d0gb(:,:)     * sqrt(sx)
      Basis%d1gb(:,:,:)   = Basis%d1gb(:,:,:)   * sqrt(sx)*sx
      Basis%d2gb(:,:,:,:) = Basis%d2gb(:,:,:,:) * sqrt(sx)*sx*sx
    ELSE
      write(out_unitp,*) ' ERROR in Scale_Basis'
      write(out_unitp,*) ' sx is too small  or ...'
      write(out_unitp,*) ' the basis is not allocated.'
      write(out_unitp,*) 'Sx', Sx
      write(out_unitp,*) 'alloc', Basis_IS_allocated(Basis)
      STOP 'ERROR in Scale_Basis'
    END IF
  END SUBROUTINE Scale_Basis

END MODULE Basis_m
