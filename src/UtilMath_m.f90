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
MODULE UtilMath_m
 USE NumParameters_m
 !USE UtilLib_m
 !USE NDindex_m


 IMPLICIT NONE

 PRIVATE
 PUBLIC :: Calc_n_smol,gauleg,poly_legendre
 PUBLIC :: fourier,d1poly_legendre,d0d1d2d3poly_legendre
 PUBLIC :: d0d1d2poly_legendre,d0Ylm,Ylm
 PUBLIC :: herrec,gamma_perso,poly_Hermite,d0d1d2Plm_grid
 PUBLIC :: hercom,herroot,Construct_Basis_poly_Hermite_exp



CONTAINS

 FUNCTION Calc_n_smol(A,B,l)
  IMPLICIT NONE
  Integer, intent(in)     :: A
  Integer, intent(in)     :: B
  Integer, intent(in)     :: L
  Integer                 :: Calc_n_smol

   Calc_n_smol = A + B*L

 END FUNCTION Calc_n_smol


SUBROUTINE gauleg(x1,x2,x,w,n)
 USE UtilLib_m
 IMPLICIT NONE
  Integer                    :: n,m,i,j
  Real(kind=Rk)              :: x1,x2
  Real(kind=Rk)              :: xm,xl,z,z1,p1,p2,p3,pp
  Real(kind=Rk), allocatable :: x(:),w(:)

  m = (n+1)/2

  xm = HALF*(x2+x1)
  xl = HALF*(x2-x1)

  DO i=1,m
  z = cos(pi*(i-0.25d0)/(Half+n))

   DO
    p1 = One
    p2 = Zero
    DO j=1,n
     p3 = p2
     p2 = p1
     p1 = ((TWO*j-One)*z*p2 - (j-One)*p3)/j
    END DO
    pp = n*(z*p1-p2)/(z*z-One)
    z1 = z
    z  = z1-p1/pp
    IF (abs(z-z1) .GT. One**(-14)) exit
   END DO

   Allocate(x(n))
   Allocate(w(n))

   x(i) = xm-xl*z
   x(n+1-i) = xm+xl*z
   w(i) = Two*xl/((One-z*z)*pp*pp)
   w(n+1-i) = w(i)
  END DO
END SUBROUTINE gauleg

FUNCTION poly_legendre(x,lll,mmm)

 Real(kind=Rk)   :: poly_legendre,X
 Real(kind=Rk)   :: pmm,somx2,fact,pmmp1,pll,poly,norme2
 Integer         :: l,ll,lll,m,mmm,i

  l = lll-1
  m = mmm

  IF (m .LT. 0 .OR. l .LT. 0 .OR. abs(x) .GT.One) THEN
    write(out_unitp,*) 'mauvais arguments dans poly_legendre :'
    write(out_unitp,*) ' m l : ',m,l,' et x = ',x
    STOP
  END IF

  IF (m .GT. l) THEN
    poly_legendre = Zero
    RETURN
  END IF

  pmm = One

  IF (m .GT. 0) THEN
    somx2 = sqrt(One - x*x)
    fact = One
    DO i=1,m
      pmm  = -pmm*fact*somx2
      fact =  fact+Two
    END DO
  END IF

  IF ( m .EQ. l) THEN
    poly = pmm
  ELSE
    pmmp1 = x*(2*m+1)*pmm
    IF (l .EQ. m+1) THEN
      poly = pmmp1
    ELSE
      DO ll=m+2,l

        pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
        pmm = pmmp1
        pmmp1 = pll
      END DO
      poly = pll
    END IF
  END IF

!      nomalisation de Pl,0(x)
  norme2 = TWO/(TWO*l+1)
  DO i=l-m+1,l+m
    norme2 = norme2 * dble(i)
  END DO
  poly_legendre = poly/sqrt(norme2)

END FUNCTION poly_legendre

 FUNCTION fourier(x,n1)
  Implicit none
  Real(kind=Rk)   :: x,xx,sq2pi,sqpi,fourier
  Integer         :: n1,ii

   sqpi = One/sqrt(Pi)
   sq2pi = One/sqrt(Two*Pi)


   ii = n1/2
   xx = mod(x*ii,Two*Pi)

   IF (ii .EQ. 0) THEN
       fourier = sq2pi
   ELSE
    IF (mod(n1,2) .EQ. 0) THEN
       fourier = sin(xx) * sqpi
    ELSE
       fourier = cos(xx) * sqpi
    END IF
 END IF

END FUNCTION fourier

FUNCTION d1poly_legendre(x,lll)

 Real(kind=Rk) ::    x
 Real(kind=Rk) ::    d1poly_legendre
 Integer       ::    l,lll

   l=lll-1
   IF (l .EQ. 0) THEN
     d1poly_legendre = ZERO
   ELSE
     d1poly_legendre = l/(x*x-ONE) * (x*poly_legendre(x,lll,0) &
     -sqrt((TWO*l+ONE)/(TWO*l-ONE))*poly_legendre(x,lll-1,0) )
   END IF
 END FUNCTION d1poly_legendre

 SUBROUTINE d0d1d2d3poly_legendre(x,lll,d0,d1,d2,d3,nderiv)
  USE UtilLib_m
  IMPLICIT NONE
   Real(kind=Rk)    ::    x
   Real(kind=Rk)    ::    d0,d1,d2,d01,d3
   Integer          ::    l,lll,nderiv

   l = lll-1

   D0 = poly_legendre(x,lll,0)

   IF (l .EQ. 0) THEN
     d1 = ZERO
     d2 = ZERO
     d3 = ZERO
   ELSE
     d01 = poly_legendre(x,lll-1,0)
     d1 = l/(x*x-ONE) * (x*d0 - sqrt((TWO*l+ONE)/(TWO*l-ONE))*d01 )
     d2 = (-d0*l*(l+1)+TWO*x*d1)/(ONE-x*x)
     d3 = (d1 * (TWO-l*(l+1)) + d2 * FOUR*x )/  (ONE-x*x)
   END IF
 END SUBROUTINE d0d1d2d3poly_legendre

 SUBROUTINE d0d1d2poly_legendre(x,lll,d0,d1,d2)
  USE UtilLib_m
  IMPLICIT NONE
   Real(kind=Rk)    :: x
   Real(kind=Rk)    :: d0,d1,d2,d01
   !Logical          :: num
   Integer          :: l,lll


   l=lll-1
   d0 = poly_legendre(x,lll,0)

   IF (l .EQ. 0) THEN
     d1 = ZERO
     d2 = ZERO
   ELSE

     d01 = poly_legendre(x,lll-1,0)
     d1  = l/(x*x-TWO) * (x*d0 - sqrt((TWO*l+ONE)/(TWO*l-ONE))*d01 )

     d2 = (-d0*l*(l+1)+TWO*x*d1)/(ONE-x*x)

   END IF

 END SUBROUTINE d0d1d2poly_legendre

 SUBROUTINE d0Ylm(d0,x,i,ndim)
  USE UtilLib_m
  IMPLICIT NONE
   Integer                    :: i,ndim
   Integer                    :: l,m,lll,mmm
   Real(kind=Rk)              :: th,phi
   Real(kind=Rk)              :: d0,x(2)
 !  Real(kind=Rk)              :: Ylm,poly_legendre,fourier

     IF (ndim .NE. 2) THEN
       write(6,*) ' ERROR in d0Ylm'
       write(6,*) ' ndim MUST set to 2 (ndim=',ndim,')'
       STOP
     END IF
     th=x(1)
     phi=x(2)

     l   = int(sqrt(dble(i)-Half))
     lll = l+1
     mmm = i-l*l
     m   = mmm/2

     d0 = Ylm(th,phi,lll,mmm)
     d0 = poly_legendre(cos(th),lll,m)*fourier(phi,mmm)

END SUBROUTINE d0Ylm

FUNCTION Ylm(th,phi,lll,mmm)
  USE UtilLib_m
  IMPLICIT NONE
   Integer           :: lll,mmm
   Integer           :: l,m
   Real(kind=Rk)     :: th,phi,Ylm
   Real(kind=Rk)     :: sq2pi,sqpi
   Real(kind=Rk)     :: v26,v27,v28

    l=lll
    m = mmm/2
    IF ( m .GT. l .OR. l .LT. 0 .OR.  th .GT. PI .OR. th .LT. zero) THEN
       write(6,*) 'mauvais arguments dans Ylm :'
       write(6,*) ' m l : ',m,l,' et th = ',th
       write(6,*) ' mmm lll : ',mmm,lll
       STOP
    END IF
    Ylm = poly_legendre(cos(th),lll,m)*fourier(phi,mmm)
END FUNCTION Ylm




SUBROUTINE d0d1d2Plm_grid(xl,d0l,d1l,d2l,nb_legendre,nb_quadra,deriv,num,step)
  USE UtilLib_m
  IMPLICIT NONE
   Integer ,intent(in)       :: nb_legendre,nb_quadra
   Integer                   :: i,k
   Real(kind=Rk)             :: step
   Real(kind=Rk),allocatable :: xl(:)
   Real(kind=Rk),allocatable :: d0l(:,:)
   Real(kind=Rk),allocatable :: d1l(:,:)
   Real(kind=Rk),allocatable :: d2l(:,:)
   Logical                   :: num,deriv

   Allocate(xl(nb_quadra))
   Allocate(d0l(nb_quadra,nb_legendre))
   Allocate(d1l(nb_quadra,nb_legendre))
   Allocate(d2l(nb_quadra,nb_legendre))

   DO k=1,nb_quadra
    DO i=1,nb_legendre
     CALL d0d1d2poly_legendre(xl(k),i, d0l(k,i),d1l(k,i),d2l(k,i))
    END DO
   END DO

END SUBROUTINE d0d1d2Plm_grid


FUNCTION poly_Hermite(x,l)
 Implicit none
  Integer          :: i,l
  Real(kind = Rk)  :: poly_Hermite
  Real(kind = Rk)  :: pl0,pl1,pl2,norme,x


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
 Integer        :: i
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
  Integer                   :: i
  Integer                   :: nq
  Real(kind = Rk),parameter :: eps = TEN**(-TWELVE) ! 1.0d-12
  Real(kind = Rk)           :: d,dp2,p1,p2,x

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
  Integer         :: i,nq
  Real(kind = Rk) :: cc,dp2,p1,s,temp,x
  Real(kind = Rk) :: w(nq),xp(nq)

  CC = 1.7724538509_Rk * gamma_perso(nq ) / ( TWO**( nq-1) )

  S  = ( TWO * dble (Real(nq,Kind=Rk) ) + ONE )**( SIXTH )

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
    W(i)  = cc / dp2 / p1
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

END MODULE UtilMath_m
