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
MODULE diago_m
!$ USE omp_lib
  IMPLICIT NONE

   PRIVATE
   PUBLIC diagonalization,Diago_Arpack

  INTERFACE diagonalization
     MODULE PROCEDURE QML_diagonalization
  END INTERFACE

   CONTAINS
!============================================================
!
!   Driver for the diagonalization
!      Default: tred2+tql2 (type_diag=2)
!            Other possibilities: Jacobi (type_diag=1) or Lapack (type_diag=3)
!            Rk: Lapack diago is not possible
!      Sort: the eigenvalues/eigenvectors:
!            sort=1:  ascending (default)
!            sort=-1: descending
!            sort=2:  ascending on the absolute eigenvalues
!     phase:
!============================================================
!
  RECURSIVE SUBROUTINE QML_diagonalization(Mat,REig,Vec,n,type_diag,sort,phase,IEig)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real64,int32
    USE NumParameters_m
    IMPLICIT NONE

    integer,       intent(in)              :: n ! when n < size(REig), only n eigenvectors and  eigenvectors are calculated
    real(kind=Rk), intent(in)              :: Mat(:,:)
    real(kind=Rk), intent(inout)           :: REig(:),Vec(:,:)
    real(kind=Rk), intent(inout), optional :: IEig(:)

    integer,          intent(in),    optional :: type_diag,sort
    logical,          intent(in),    optional :: phase


    !local variables
    integer            :: type_diag_loc
    integer            :: type_diag_default = 2 ! tred+tql
    !                                    Jacobi tred+tql DSYEV  DGEEV Lanczos
    integer, parameter :: list_type(7) = [1,    2,       3,377, 4,477,   5]

    real(kind=Rk), allocatable :: trav(:),Mat_save(:,:)
    integer :: n_size,n_vect

    !for lapack
    integer              :: i
    integer              :: lwork ,lda ,ldvr ,ierr
    integer(kind=int32)  :: n4,lwork4,lda4,ldvr4,ierr4
    real(kind=Rk), allocatable :: work(:)
    real(kind=Rk), allocatable :: IEig_loc(:)

    real(kind=Rk) :: dummy(1,1)





    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='QML_diagonalization'
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    !-----------------------------------------------------------

    n_size = size(REig)
    IF (n_size /= size(Mat,dim=1) .OR. n_size /= size(Vec,dim=1)) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The matrix or vector sizes are not consistant'
      write(out_unitp,*) '   size(REig):     ',size(REig)
      write(out_unitp,*) '   size(Mat):      ',size(Mat,dim=1)
      write(out_unitp,*) '   size(Vec):      ',size(Vec,dim=1)
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in QML_diagonalization: The matrix or vector sizes are not consistant.'
    END IF
    IF (n < 1) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,"(a,i0,a)") ' n < 1. It MUST be in the range [1,',n_size,']'
      write(out_unitp,*) '   n:              ',n
      write(out_unitp,*) '   size(REig):     ',size(REig)
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in QML_diagonalization: The matrix or vector sizes are not consistant.'
    END IF
    n_vect = min(n,n_size)

    IF (present(type_diag)) THEN
      type_diag_loc = type_diag
    ELSE
      type_diag_loc = type_diag_default
    END IF

    !when lapack is used and Rk /= real64 (not a double)
    IF (Rk /= real64 .AND. type_diag_loc == 3) type_diag_loc = type_diag_default

!#if __LAPACK != 1
    IF (count([3,377,395] == type_diag_loc) == 1) type_diag_loc = type_diag_default
    IF (count([4,477] == type_diag_loc) == 1) THEN
      !type_diag_loc = 0
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The diagonalization of non-symmetric needs LAPACK.'
      write(out_unitp,*) '  Try to link LAPACK with the code (use LAPACK=1 in the makfile).'
      write(out_unitp,*) '   type_diag:      ',type_diag_loc
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in QML_diagonalization: Problem with non-symmetric matrix.'
    END IF
!#endif

    IF (count(list_type == type_diag_loc) == 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' type_diag is out-of-range.'
      write(out_unitp,*) '   type_diag:      ',type_diag_loc
      write(out_unitp,*) '   Possible values:',list_type(:)
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in QML_diagonalization: type_diag is out-of-range.'
    END IF


    SELECT CASE (type_diag_loc)
    CASE(1) ! jacobi
      IF (debug) write(out_unitp,*) 'Jacobi (symmetric)'
      allocate(Mat_save(n,n))
      Mat_save = Mat ! save mat

      CALL QML_JACOBI2(Mat_save,n,REig,Vec)

      deallocate(Mat_save)
    CASE (2) ! tred+tql
      IF (debug) write(out_unitp,*) 'tred+tql, new version (symmetric)'
      allocate(trav(n))

      Vec = Mat
      CALL QML_TRED2_EISPACK(Vec,n,n,REig,trav)
      CALL QML_TQLI_EISPACK(REig,trav,n,n,Vec)

      deallocate(trav)
    CASE(3,377) ! lapack77
      IF (debug) write(out_unitp,*) 'lapack77: DSYEV (symmetric)'


      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makefile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in QML_diagonalization: LAPACK case impossible'
!#endif

    CASE(4,477) ! lapack77 (non-symmetric)

      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) '  LAPACK is not linked (LAPACK=0 in the makefile).'
      write(out_unitp,*) '  The program should not reach the LAPACK case.'
      write(out_unitp,*) '  => Probabely, wrong type_diag_default.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in QML_diagonalization: LAPACK case impossible'

    CASE(5) ! lanczos

      CALL QML_Lanczos(Mat,n_vect,REig,Vec,epsi=ONETENTH**6,max_it=100)

    CASE DEFAULT
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' The default CASE is not defined.'
      write(out_unitp,*) '  => CHECK the fortran!!'
      STOP 'ERROR in QML_diagonalization: default case impossible'
    END SELECT


    IF (present(sort)) THEN
        SELECT CASE (sort)
        CASE(1)
          CALL QML_sort(REig,Vec)
          CALL QML_rota_denerated(REig,Vec)
        CASE(-1)
          REig = -REig
          CALL QML_sort(REig,Vec)
          REig = -REig
          CALL QML_rota_denerated(REig,Vec)
        CASE(2)
          CALL QML_sort_abs(REig,Vec)
        CASE DEFAULT ! no sort
          CONTINUE
        END SELECT
    ELSE
      CALL QML_sort(REig,Vec)
      CALL QML_rota_denerated(REig,Vec)
    END IF

    IF (present(phase)) THEN
      IF (phase) CALL QML_Unique_phase(Vec)
    ELSE
      CALL QML_Unique_phase(Vec)
    END IF

  END SUBROUTINE QML_diagonalization

  SUBROUTINE QML_JACOBI2(A,N,D,V)
      USE NumParameters_m
      IMPLICIT NONE

      integer            :: N
      real(kind=Rk)   :: A(N,N),V(N,N),D(N)


      integer, parameter :: max_it = 500
      real(kind=Rk)   :: B(N),Z(N)


      real(kind=Rk)   :: h,t,g,sm,tresh,tau,s,theta,c

      integer            :: i,j,iq,nrot,ip

!     V(:,:) = Id(:,:)
      V(:,:) = ZERO
      DO IP=1,N
        V(IP,IP)=ONE
      END DO

!     initialization
      DO IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=ZERO
      END DO

      NROT=0
      DO I=1,max_it ! main loop

        ! SM value
        SM = ZERO
        DO IP=1,N-1
          DO IQ=IP+1,N
            SM = SM+abs(A(IP,IQ))
          END DO
        END DO
        IF(SM == ZERO)RETURN

        ! TRESH value
        IF(I < 4)THEN
          TRESH = TWOTENTHS*SM/N**2
        ELSE
          TRESH = ZERO
        ENDIF

        DO IP=1,N-1
          DO IQ=IP+1,N
            G = HUNDRED*abs(A(IP,IQ))
            IF ( I > 4 .AND. ABS(D(IP))+G == ABS(D(IP))                 &
               .AND. ABS(D(IQ))+G == ABS(D(IQ)) ) THEN
              A(IP,IQ)=ZERO
            ELSE IF ( ABS(A(IP,IQ)) > TRESH ) THEN
              H=D(IQ)-D(IP)
              IF ( ABS(H)+G == ABS(H) ) THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=HALF*H/A(IP,IQ)
                T=ONE/(ABS(THETA)+sqrt(ONE+THETA**2))
                IF ( THETA < ZERO) T=-T
              ENDIF
              C=ONE/sqrt(ONE+T**2)
              S=T*C
              TAU=S/(ONE+C)

              H=T*A(IP,IQ)

              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=ZERO
              DO J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
              END DO
              DO J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
              END DO
              DO J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
              END DO
              DO J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
              END DO
              NROT=NROT+1
            ENDIF
          END DO
        END DO

        DO IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=ZERO
        END DO

      END DO ! end main loop

      write(out_unitp,*) max_it,' iterations should never happen'
      STOP

  end subroutine QML_JACOBI2
!
!============================================================
!
!   diagonalisation trigonalisation puis diagonalisation
!
!============================================================
!
  SUBROUTINE QML_TRED2_EISPACK(A,N,NP,D,E)
      USE NumParameters_m
      IMPLICIT NONE

      integer          :: N,NP
      real(kind=Rk) :: A(NP,NP),D(NP),E(NP)

      !local variables
      integer          :: I,J,K,L
      real(kind=Rk) :: F,G,H,HH,SCALE

      IF(N.GT.1)THEN
        DO 18 I=N,2,-1
          L=I-1
          H=ZERO
          SCALE=ZERO
          IF(L.GT.1)THEN
            DO 11 K=1,L
              SCALE=SCALE+ABS(A(I,K))
11          CONTINUE
            IF(SCALE .EQ. ZERO)THEN
              E(I)=A(I,L)
            ELSE
              DO 12 K=1,L
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
12            CONTINUE
              F=A(I,L)
              G=-SIGN(SQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=ZERO
              DO 15 J=1,L
                A(J,I)=A(I,J)/H
                G=ZERO
                DO 13 K=1,J
                  G=G+A(J,K)*A(I,K)
13              CONTINUE
                IF(L.GT.J)THEN
                  DO 14 K=J+1,L
                    G=G+A(K,J)*A(I,K)
14                CONTINUE
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
15            CONTINUE
              HH=F/(H+H)
              DO 17 J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO 16 K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
16              CONTINUE
17            CONTINUE
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
18      CONTINUE
      ENDIF
      D(1)=ZERO
      E(1)=ZERO
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE. ZERO)THEN
          DO 21 J=1,L
            G=ZERO
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=ONE
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=ZERO
            A(J,I)=ZERO
22        CONTINUE
        ENDIF
23    CONTINUE
      RETURN
  END SUBROUTINE QML_TRED2_EISPACK

  SUBROUTINE QML_TQLI_EISPACK(D,E,N,NP,Z)
      USE NumParameters_m
      IMPLICIT NONE

      integer          :: N,NP
      real(kind=Rk) :: D(NP),E(NP),Z(NP,NP)

      !local variables
      integer          :: I,K,L,M,ITER
      real(kind=Rk) :: G,R,S,C,P,F,B,DD

      IF (N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
11      CONTINUE
        E(N)=ZERO
        DO 15 L=1,N
          ITER=0
1         DO 12 M=L,N-1
            DD=ABS(D(M))+ABS(D(M+1))
            IF (ABS(E(M))+DD.EQ.DD) GO TO 2
12        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30) STOP 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(TWO*E(L))
            R=SQRT(G**2+ONE)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=ONE
            C=ONE
            P=ZERO
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(ABS(F).GE.ABS(G))THEN
                C=G/F
                R=SQRT(C**2+ONE)
                E(I+1)=F*R
                S=ONE/R
                C=C*S
              ELSE
                S=F/G
                R=SQRT(S**2+ONE)
                E(I+1)=G*R
                C=ONE/R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+TWO*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
13            CONTINUE
14          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=ZERO
            GO TO 1
          ENDIF
15      CONTINUE
      ENDIF
      RETURN
  END SUBROUTINE QML_TQLI_EISPACK

!
!============================================================
!
!   Sort the eigenvalues and the eigenvectors
!
!============================================================
!
  SUBROUTINE QML_sort_tab(nb_niv,ene,max_niv)
      USE NumParameters_m
      IMPLICIT NONE

      integer nb_niv,max_niv
      real(kind=Rk) ene(max_niv)
      real(kind=Rk) a

        integer i,j,k

      DO i=1,nb_niv
      DO j=i+1,nb_niv
       IF (ene(i) .GT. ene(j)) THEN
!             permutation
          a=ene(i)
          ene(i)=ene(j)
          ene(j)=a
        END IF
      END DO
      END DO

  end subroutine QML_sort_tab

  SUBROUTINE QML_sort(ene,psi)
      USE NumParameters_m
      IMPLICIT NONE

      real(kind=Rk), intent(inout) :: ene(:),psi(:,:)

      real(kind=Rk) :: a
      integer          :: i,j,k

      DO i=1,size(ene)
      DO j=i+1,size(ene)
       IF (ene(i) > ene(j)) THEN
!	      permutation
          a=ene(i)
          ene(i)=ene(j)
          ene(j)=a
          DO k=1,size(psi,dim=1)
            a=psi(k,i)
            psi(k,i)=psi(k,j)
            psi(k,j)=a
          END DO
        END IF
      END DO
      END DO

  END SUBROUTINE QML_sort

  SUBROUTINE QML_sort_abs(ene,psi)
      USE NumParameters_m
      IMPLICIT NONE

      real(kind=Rk), intent(inout) :: ene(:),psi(:,:)

      real(kind=Rk) :: a
      integer          :: i,j,k


        DO i=1,size(ene)
          DO j=i+1,size(ene)
            IF (abs(ene(i)) > abs(ene(j))) THEN
!             permutation
              a=ene(i)
              ene(i)=ene(j)
              ene(j)=a
              DO k=1,size(psi,dim=1)
                a=psi(k,i)
                psi(k,i)=psi(k,j)
                psi(k,j)=a
              END DO
            END IF
          END DO
        END DO

  end subroutine QML_sort_abs
!
!============================================================
!
!   Change the phase of Vec(:,i) shuch its largest coeficient is positive
!
!============================================================
!
  SUBROUTINE QML_Unique_phase(Vec)
      USE NumParameters_m
      IMPLICIT NONE

      real(kind=Rk), intent(inout) :: Vec(:,:)

      integer          :: i,jloc

      DO i=lbound(Vec,dim=2),ubound(Vec,dim=2)
        jloc           = maxloc(abs(Vec(:,i)),dim=1)
        IF (abs(Vec(jloc,i)) < ONETENTH**6 ) CYCLE
        IF (Vec(jloc,i) < ZERO) Vec(:,i) = -Vec(:,i)
      END DO

      END SUBROUTINE QML_Unique_phase
!=====================================================================
!
!   c_new(:,i) =  cos(th) c(:,i) + sin(th) c(:,j)
!   c_new(:,j) = -sin(th) c(:,j) + cos(th) c(:,j)
!
!    the angle is obtained such ...
!
!      it is working only if 2 vectors are degenerated !!!!
!
!=====================================================================
  SUBROUTINE QML_rota_denerated(v,c)
      USE NumParameters_m
      IMPLICIT NONE

      real (kind=Rk), intent(in)    :: v(:)
      real (kind=Rk), intent(inout) :: c(:,:)

      integer           :: i,j,k,kloc
      real (kind=Rk) :: ai,aj,norm,cc,ss

      real (kind=Rk), parameter :: epsi = ONETENTH**10

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING QML_rota_denerated'
      write(out_unitp,*) 'v',v
      !write(out_unitp,*) 'c',c
      END IF
!---------------------------------------------------------------------
      DO i=1,size(v)-1

        j = i+1
        IF ( abs(v(i)-v(j)) < epsi) THEN
          !write(6,*) 'i,j',i,j
          !write(6,*) 'vec i',c(:,i)
          !write(6,*) 'vec j',c(:,j)


          kloc = maxloc((c(:,i)**2+c(:,j)**2),dim=1)

          cc   =  c(kloc,i)
          ss   = -c(kloc,j)
          !write(6,*) i,j,'cos sin',kloc,cc,ss
          norm = sqrt(cc*cc+ss*ss)
          cc   = cc/norm
          ss   = ss/norm
          !write(6,*) i,j,'cos sin',cc,ss
          IF (abs(cc) < epsi .OR. abs(ss) < epsi) CYCLE

          DO k=1,size(c,dim=1)
           ai = c(k,i)
           aj = c(k,j)

           c(k,i) =  cc * ai + ss * aj
           c(k,j) = -ss * ai + cc * aj

          END DO

        END IF
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'new c',c
      write(out_unitp,*) 'END QML_rota_denerated'
      END IF
!---------------------------------------------------------------------

  end subroutine QML_rota_denerated

!============================================================
!
!      *******       Extension to complex symmetrical matrices of
!      *******       the <Tql2> Eispack routine, implemented by
!      *******       Claude Leforestier.
!
!
!     on input-
!        n is the order of the matrix,
!        d contains the diagonal elements of the input matrix,
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e(1) is arbitrary,
!        z should have been initialized to the identity matrix.
!
!      on output-
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1,
!        e has been destroyed,
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues,
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after nbiter iterations.
!
!============================================================
!
  Subroutine QML_cTql2(nZ,n,D,E,Z,ierr)
      USE NumParameters_m
      IMPLICIT NONE

      Integer :: i,j,k,l,m,n,ii,Nbiter,nm,mml,nZ,ierr
      Data Nbiter/60/
      complex(kind=Rk) :: D(n),E(n),Z(nZ,n)
      real (kind=Rk) :: machep,norm1,norm2,rsign
      complex(kind=Rk) :: b,c,f,g,p,r,s
!     complex(kind=Rk) :: b,c,f,g,p,r,s,zero,one,two
!     Data zero/(0.,0.)/,one/(1.,0.)/,two/(2.,0.)/,Nbiter/60/


      machep=10.d-16
      ierr=0
!     initialize z to e one.
      do i=2,n
         e(i-1)=e(i)
      enddo
      e(n)=zero
      do 240 l=1,n
      j=0
!     ********** look for small sub-diagonal element **********
  105 do 110 m=l,n
      if(m.eq.n) goto 120
      norm1=abs(real(e(m),kind=Rk))+abs(aimag(e(m)))
      norm2=abs(real(d(m),kind=Rk))+abs(aimag(d(m))) +               &
            abs(real(d(m+1),kind=Rk))+abs(aimag(d(m+1)))
      if(norm1.le.machep*norm2) goto 120
  110 continue
  120 p=d(l)
      if(m.eq.l) goto 240
      if(j.eq.nbiter) goto 1000
      j=j+1
!     ********** form shift **********
      g=(d(l+1)-p)/(TWO*e(l))
      r=sqrt(g*g+ONE)
      rsign=1.
      if(real(g,kind=Rk).lt.0.) rsign=-ONE
      g=d(m)-p+e(l)/(g+rsign*r)
      s=ONE
      c=ONE
      p=ZERO
      mml=m-l
!     ********** for i=m-1 step -1 until l do -- **********
      do 200 ii=1,mml
      i=m-ii
      f=s*e(i)
      b=c*e(i)
      norm1=f*conjg(f)
      norm2=g*conjg(g)
      if(norm1.lt.norm2) goto 150
      c=g/f
      r=sqrt(c*c+ONE)
      e(i+1)=f*r
      s=ONE/r
      c=c*s
      go to 160
  150 s=f/g
      r=sqrt(s*s+ONE)
      e(i+1)=g*r
      c=ONE/r
      s=s*c
  160 g=d(i+1)-p
      r=(d(i)-g)*s+TWO*c*b
      p=s*r
      d(i+1)=g+p
      g=c*r-b
!     ********** form vector **********
      do k=1,n
      f=z(k,i+1)
      z(k,i+1)=s*z(k,i)+c*f
      z(k,i)=c*z(k,i)-s*f
      enddo
  200 continue
      d(l)=d(l)-p
      e(l)=g
      e(m)=ZERO
      go to 105
  240 continue
!     ********** order eigenvalues and eigenvectors **********
      do 300 ii=2,n
      i=ii-1
      k=i
      p=d(i)
      do 260 j=ii,n
      if(real(d(j),kind=Rk).ge.real(p,kind=Rk)) goto 260
      k=j
      p=d(j)
  260 continue
      if(k.eq.i) goto 300
      d(k)=d(i)
      d(i)=p
      do m=1,n
         p=z(m,i)
         z(m,i)=z(m,k)
         z(m,k)=p
      enddo
  300 continue
      goto 1001
!     ********** set error -- no convergence to an
!                eigenvalue after Nbiter iterations **********
 1000 write (out_unitp,1010) l
 1010 format(//10x,'$$$ <Cmtql2> return code :',i5,' $$$')
      stop
 1001 Return
  end subroutine QML_cTql2

  Subroutine QML_cTred2(nm,n,A,d,e,Z)
      USE NumParameters_m
      IMPLICIT NONE

      INTEGER I,II,J,JP1,K,L,N,NM
      complex(kind=Rk) A(nm,n),Z(nm,n),d(n),e(n)                     &
                ,f,g,h,hh
      real (kind=Rk) :: scale,rsign
!     real (kind=Rk) :: one,scale,rsign,zero
!     Data zero/0./,one/1./

      do i=1,n
         do j=1,i
            Z(i,j)=A(i,j)
         enddo
      enddo
!
!     FOR I=N STEP -1 UNTIL 2 DO --
      do ii=2,n
         I=N+2-II
         L=I-1
         H=ZERO
         SCALE=ZERO
         IF(L.LT.2) GOTO 130
!     SCALE ROW (ALGOL TOL THEN NOT NEEDED)
         DO 120 K=1,L
  120    SCALE=SCALE+ABS(Z(I,K))
         IF(SCALE.NE.ZERO) GOTO 140
  130    E(I)=Z(I,L)
         GOTO 290
  140    DO 150 K=1,L
         Z(I,K)=Z(I,K)/SCALE
         H=H+Z(I,K)*Z(I,K)
  150 CONTINUE
         F=Z(I,L)
         rsign=ONE
         if( real(F,kind=Rk).lt.0.) rsign=-ONE
         G=-rsign*SQRT(H)
         E(I)=SCALE*G
         H=H-F*G
         Z(I,L)=F-G
         F=ZERO
         DO 240 J=1,L
         Z(J,I)=Z(I,J)/(SCALE*H)
         G=ZERO
!     FORM ELEMENT OF A*U
         DO 180 K=1,J
  180    G=G+Z(J,K)*Z(I,K)
         JP1=J+1
         IF(L.LT.JP1) GOTO 220
         DO 200 K=JP1,L
  200    G=G+Z(K,J)*Z(I,K)
!     FORM ELEMENT OF PP
  220    E(J)=G/H
         F=F+E(J)*Z(I,J)
  240    CONTINUE
         HH=F/(H+H)
!     FORM REDUCED A
         DO 260 J=1,L
         F=Z(I,J)
         G=E(J)-HH*F
         E(J)=G
         DO 260 K=1,J
         Z(J,K)=Z(J,K)-F*E(K)-G*Z(I,K)
  260 CONTINUE
         DO 280 K=1,L
  280    Z(I,K)=SCALE*Z(I,K)
  290    D(I) = H
      enddo
  320 D(1)=ZERO
      E(1)=ZERO
!     ACCUMULATION OF TRANSFORMATION MATRICES
      DO 500 I=1,N
      L=I-1
      IF(abs(D(I)).EQ.ZERO) GOTO 380
      DO 360 J=1,L
      G=ZERO
      DO 340 K=1,L
  340 G=G+Z(I,K)*Z(K,J)
      DO 360 K=1,L
      Z(K,J)=Z(K,J)-G*Z(K,I)
  360 CONTINUE
  380 D(I)=Z(I,I)
      Z(I,I)=ONE
      IF(L.LT.1) GOTO 500
      DO 400 J=1,L
      Z(I,J)=ZERO
      Z(J,I)=ZERO
  400 CONTINUE
  500 CONTINUE
      RETURN
  end subroutine QML_cTred2

  SUBROUTINE QML_Lanczos(A,n_vect,D,V,epsi,max_it)
      USE NumParameters_m
      USE UtilLib_m
      IMPLICIT NONE

      integer,          intent(in)    :: n_vect
      real(kind=Rk), intent(in)    :: A(:,:)
      real(kind=Rk), intent(inout) :: V(:,:),D(:)
      real(kind=Rk), intent(in)    :: epsi
      integer,          intent(in)    :: max_it



      real(kind=Rk), allocatable :: Krylov_vectors(:,:)
      real(kind=Rk), allocatable :: M_Krylov(:,:)
      real(kind=Rk), allocatable :: V_Krylov(:,:)
      real(kind=Rk), allocatable :: E0_Krylov(:)
      real(kind=Rk), allocatable :: E1_Krylov(:)
      real(kind=Rk)              :: y,maxdiff

      integer :: i,it,n_size

      ! Begin Lanczos scheme
      n_size = size(A,dim=1)
      write(6,*) 'shape(A)',shape(A)
      write(6,*) 'shape(V)',shape(V)
      write(6,*) 'shape(D)',shape(D)
      write(6,*) 'n_size',n_size
      write(6,*) 'n_vect',n_vect
      write(6,*) 'max_it',max_it
      write(6,*) 'epsi',epsi

      allocate(Krylov_vectors(n_size,0:max_it))
      Krylov_vectors(:,:) = ZERO
      CALL random_number(Krylov_vectors(1:n_vect,0))
      Krylov_vectors(:,0) = Krylov_vectors(:,0) /                               &
                      sqrt(dot_product(Krylov_vectors(:,0),Krylov_vectors(:,0)))
      write(6,*) 'Krylov_vectors (guess)',Krylov_vectors(:,0)


      allocate(M_Krylov(max_it,max_it))
      allocate(V_Krylov(max_it,max_it))
      allocate(E0_Krylov(max_it))
      allocate(E1_Krylov(max_it))
      maxdiff = huge(ONE)

      DO it=1,max_it

         Krylov_vectors(:,it) = matmul(A,Krylov_vectors(:,it-1))

         DO i=0,it-1
            M_Krylov(i+1, it) = dot_product(Krylov_vectors(:,i),Krylov_vectors(:,it))
            M_Krylov(it, i+1) = M_Krylov(i+1,it)
         END DO
         CALL Write_RMat(M_Krylov(1:it,1:it),out_unitp,5)

         ! Orthogonalize vectors
         DO i=0,it-1
            y = dot_product(Krylov_vectors(:,i),Krylov_vectors(:,it))
            Krylov_vectors(:,it) = Krylov_vectors(:,it) - y * Krylov_vectors(:,i)
         END DO
         ! Normalize vector
          Krylov_vectors(:,it) =   Krylov_vectors(:,it) /                      &
                  sqrt(dot_product(Krylov_vectors(:,it),Krylov_vectors(:,it)))

         IF (it >= n_vect) THEN
            CALL diagonalization(M_Krylov(1:it,1:it),E1_Krylov(1:it),V_Krylov(1:it,1:it),it,3,1,.FALSE.)
            write(6,*) 'it eig',it,E1_Krylov(1:n_vect)
         END IF
         IF (it > n_vect) THEN
            maxdiff = maxval(abs(E1_Krylov(1:n_vect) - E0_Krylov(1:n_vect)))
            E0_Krylov(1:n_vect) = E1_Krylov(1:n_vect)
         ELSE
            maxdiff = huge(ONE)
         END IF
         write(6,*) 'it maxdiff,epsi,exit?',it,maxdiff,epsi,(maxdiff < epsi)

         IF (maxdiff < epsi) EXIT

      END DO
      stop
  end subroutine QML_Lanczos


  SUBROUTINE Diago_Arpack(psi,Ene,nb_diago,max_diago,max_it,para_H)
    USE NumParameters_m

    IMPLICIT NONE

    !----- Operator: Hamiltonian ----------------------------
    real (kind=Rk)            :: para_H(:,:)
    integer                   :: nb_diago,max_diago,max_it
    real (kind=Rk)            :: psi(:,:)
    real (kind=Rk)            :: Ene(max_diago)
    logical, parameter        :: temp_logi = .false.
    !------ working parameters --------------------------------
    real (kind=Rk)              :: ZPE
    real (kind=Rk), allocatable :: Evec(:)

!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%
    integer  (kind=4)            :: maxn, maxnev, maxncv, ldv
!     %--------------%
!     | Local Arrays |
!     %--------------%
    real(kind=Rk), allocatable    :: v(:,:),workl(:),workd(:),workev(:),d(:,:),resid(:),ax(:)
    logical, allocatable          :: select(:)
    integer  (kind=4)             :: iparam(11), ipntr(14)

!     %---------------%
!     | Local Scalars |
!     %---------------%
    character (len=1)             :: bmat
    character (len=2)             :: which
    integer (kind=4)              :: ido, n, nev, ncv, lworkl, info, ierr, j
    integer (kind=4)              :: nconv, maxitr, mode, ishfts
    logical                       :: rvec,first
    logical                       :: if_deq0=.FALSE. ! for MPI
    real(kind=Rk)                 :: tol, sigmar, sigmai
    integer (kind=4)              :: it
!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%
    real(kind=Rk)                 ::  dnrm2,dlapy2
    external                      :: dnrm2, dlapy2, daxpy

!----- for debuging --------------------------------------------------
    integer                       :: err_mem,memory
    character (len=*), parameter  :: name_sub='Diago_Arpack'
    !logical, parameter           :: debug=.FALSE.
    logical, parameter            :: debug=.TRUE.
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) ' nb_diago',nb_diago
      write(out_unitp,*) ' max_diago',max_diago
      write(out_unitp,*)
      !  CALL flush_perso(out_unitp)
    END IF
!-----------------------------------------------------------

!------ initialization -------------------------------------

    n = size(Para_H, dim=1)

    IF (nb_diago == 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' nb_diago=0 is not possible with ARPACK'
      STOP
    END IF

    nev =  int(nb_diago,kind=4)
    ncv =  2*nev+10
    ncv =  3*nev+10
      !was ncv =  2*nev+10
    ncv = MIN(N,ncv)  ! prevent infor=-3 error:
                      ! NCV must be greater than NEV and less than or equal to N.
    maxn   = n
    ldv    = maxn
    maxnev = nev
    maxncv = ncv

    allocate(select(maxncv))
    allocate(ax(maxn))
    allocate(d(maxncv,3))
    allocate(resid(maxn)) !< initial guess vector when infor=1
    allocate(v(ldv,maxncv))
    allocate(workd(3*maxn))
    allocate(workev(3*maxncv))
    allocate(workl(3*maxncv*maxncv+6*maxncv))
    bmat  = 'I'
    which = 'SR'

!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%

    lworkl = 3*ncv**2+6*ncv
    tol    = ZERO
    !tol    = ONETENTH**6

    ido    = 0
    info   = 0

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%

    ishfts = 1
    maxitr = int(min(max_diago,max_it),kind=4)
    mode   = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
    it = 0
    DO
      it = it +1
!       %---------------------------------------------%
!       | Repeatedly call the routine DNAUPD and take |
!       | actions indicated by parameter IDO until    |
!       | either convergence is indicated or maxitr   |
!       | has been exceeded.                          |
!       %---------------------------------------------%

      write(out_unitp,*) 'it,ido,n,tol,ncv,iparam,info,lworkl',it,ido,n,tol,ncv,info,lworkl
      write(out_unitp,*) "shape(iparam),=",shape(iparam)
      write(out_unitp,*) "shape(workd,)=,",shape(workd)
      write(out_unitp,*) "shape(workl),=",shape(workl)
      write(out_unitp,*) "shape(v),=",shape(v)

      call dnaupd(ido, bmat, n, which, nev, tol, resid,                 &
                ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, info)

      write(out_unitp,*) 'test'
      write(out_unitp,*) 'it,ido,n,tol,ncv,iparam,info,lworkl',it,ido,n,tol,ncv,info,lworkl
      write(out_unitp,*) "shape(iparam),=",shape(iparam)
      write(out_unitp,*) "shape(workd,)=,",shape(workd)
      write(out_unitp,*) "shape(workl),=",shape(workl)
      write(out_unitp,*) "shape(v),=",shape(v)
      write(out_unitp,*) 'it,ido',it,ido
      write(out_unitp,*) ",shape(workd(ipntr(2):ipntr(2)-1+n)),=",shape(workd(ipntr(2):ipntr(2)-1+n))
      write(out_unitp,*) "shape(para_H),=",shape(para_H)
      Stop "RObert"
      IF (abs(ido) /= 1) EXIT
!       %--------------------------------------%
!       | Perform matrix vector multiplication |
!       |              y <--- OP*x             |
!       | The user should supply his/her own   |
!       | matrix vector multiplication routine |
!       | here that takes workd(ipntr(1)) as   |
!       | the input, and return the result to  |
!       | workd(ipntr(2)).                     |
!       %--------------------------------------%
!a rajouter H:psi
!workd(ipntr(2):ipntr(2)-1+n) = H:workd(ipntr(1):ipntr(1)-1+n)



     workd(ipntr(2):ipntr(2)-1+n) = Matmul(para_H,workd(ipntr(2):ipntr(2)-1+n))

     write(out_unitp,*) 'Arpack <psi H psi>:',                                        &
                  dot_product(workd(ipntr(1):ipntr(1)-1+n),workd(ipntr(2):ipntr(2)-1+n))

  END DO


!---------------------------------------------------------------------------------------

!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

  IF ( info .lt. 0 ) THEN

!       %--------------------------%
!       | Error message. Check the |
!       | documentation in DSAUPD. |
!       %--------------------------%

    write(out_unitp,*)
    write(out_unitp,*) ' Error with _saupd, info = ', info
    write(out_unitp,*) ' Check documentation in _naupd '
    write(out_unitp,*) ' '
        !write(out_unitp,*) ' Error with _saupd, info = ', info

  ELSE

!       %-------------------------------------------%
!       | No fatal errors occurred.                 |
!       | Post-Process using DSEUPD.                |
!       |                                           |
!       | Computed eigenvalues may be extracted.    |
!       |                                           |
!       | Eigenvectors may also be computed now if  |
!       | desired.  (indicated by rvec = .true.)    |
!       %-------------------------------------------%

    rvec = .true.

    call dneupd(rvec, 'A', select, d, d(1,2), v, ldv,                           &
                      sigmar, sigmai, workev, bmat, n, which, nev, tol,    &
               resid, ncv, v, ldv, iparam, ipntr, workd, workl,lworkl, ierr )

!       %-----------------------------------------------%
!       | The real part of the eigenvalue is returned   |
!       | in the first column of the two dimensional    |
!       | array D, and the imaginary part is returned   |
!       | in the second column of D.  The corresponding |
!       | eigenvectors are returned in the first NEV    |
!       | columns of the two dimensional array V if     |
!       | requested.  Otherwise, an orthogonal basis    |
!       | for the invariant subspace corresponding to   |
!       | the eigenvalues in D is returned in V.        |
!       %-----------------------------------------------%
    IF ( ierr .ne. 0) THEN

!         %------------------------------------%
!         | Error condition:                   |
!         | Check the documentation of _neupd  |
!         %------------------------------------%

      write(out_unitp,*)
      write(out_unitp,*) ' Error with _seupd, info = ', ierr
      write(out_unitp,*) ' Check the documentation of _neupd. '
      write(out_unitp,*)

    ELSE

      first  = .true.
      nconv  = iparam(5)

      DO j=1, nconv
!           %---------------------------%
!           | Compute the residual norm |
!           |                           |
!           |   ||  A*x - lambda*x ||   |
!           |                           |
!           | for the NCONV accurately  |
!           | computed eigenvalues and  |
!           | eigenvectors.  (iparam(5) |
!           | indicates how many are    |
!           | accurate to the requested |
!           | tolerance)                |
!           %---------------------------%

        IF  (temp_logi) THEN

!             %--------------------%
!             | Ritz value is real |
!             %--------------------%
        Ax = Matmul(para_H,v(:,j))
              !CALL sub_OpV1_TO_V2_Arpack(v(:,j),ax,psi(j),Hpsi_loc,                    &
                                  ! para_H,cplxE,para_propa,int(n,4))
        call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
        d(j,3) = dnrm2(n, ax, 1)
        d(j,3) = d(j,3) / abs(d(j,1))
        Ene(j)          = d(j,1)

      ELSE IF (first) THEN

!             %------------------------%
!             | Ritz value is complex. |
!             | Residual of one Ritz   |
!             | value of the conjugate |
!             | pair is computed.      |
!             %------------------------%
              !call av(nx, v(1,j), ax)
         Ax= matmul(para_H,v(:,j))
              !CALL sub_OpV1_TO_V2_Arpack(v(:,j),ax,psi(j),Hpsi_loc,                    &
                                  !       para_H,cplxE,para_propa,int(n,4))

          call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
          call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
          d(j,3) = dnrm2(n, ax, 1)
          Ene(j) = d(j,1)

          Ax = matmul(para_H,v(:,j+1))
        !CALL sub_OpV1_TO_V2_Arpack(v(:,j+1),ax,psi(j+1),Hpsi_loc,&
                                      !  para_H,cplxE,para_propa,int(n,4))

          Ene(j+1)  = d(j,1)
          call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
          call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
          d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
          d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
          d(j+1,3) = d(j,3)
          first = .false.
      ELSE
              first = .true.
    END IF

    IF (debug) write(out_unitp,*) 'j,cplx ene ?',j,Ene(j)

    END DO ! for j=1, nconv

!         %-----------------------------%
!         | Display computed residuals. |
!         %-----------------------------%

  END IF ! for ierr .ne. 0

!       %-------------------------------------------%
!       | Print additional convergence information. |
!       %-------------------------------------------%

  IF ( info .eq. 1) THEN
    print *, ' '
    print *, ' Maximum number of iterations reached.'
    print *, ' '
  ELSE IF ( info .eq. 3) THEN
    print *, ' '
    print *, ' No shifts could be applied during implicit',    &
                   ' Arnoldi update, try increasing NCV.'
    print *, ' '
  END IF

  print *, ' '
  print *, ' _NDRV1 '
  print *, ' ====== '
  print *, ' '
  print *, ' Size of the matrix is ', n
  print *, ' The number of Ritz values requested is ', nev
  print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
  print *, ' What portion of the spectrum: ', which
  print *, ' The number of converged Ritz values is ',nconv
  print *, ' The number of Implicit Arnoldi update',             &
                 ' iterations taken is ', iparam(3)
  print *, ' The number of OP*x is ', iparam(9)
  print *, ' The convergence criterion is ', tol
  print *, ' '

 END IF ! for info .lt. 0

 nb_diago = nconv

      !CALL trie_psi(psi,Ene,nb_diago)

 IF(allocated(v)) deallocate(v)
 deallocate(workl)
 deallocate(workd)
 deallocate(d)
 deallocate(resid)
 deallocate(ax)
 deallocate(select)

!     %---------------------------%
!     | Done with program dndrv1. |
!     %---------------------------%

!----------------------------------------------------------
 IF (debug) THEN
   ZPE = minval(Evec)
   DO j=1,nb_diago
     write(out_unitp,*) j,Ene(j),(Evec(j)-ZPE)
   END DO
   write(out_unitp,*) 'END ',name_sub
 END IF
!----------------------------------------------------------
END SUBROUTINE Diago_Arpack


END MODULE diago_m
