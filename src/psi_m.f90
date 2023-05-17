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
module psi_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  implicit none

  private
  TYPE :: Vec_smol_t
    Real(kind=Rk), allocatable    :: Smol_G_l(:)
    Real(kind=Rk), allocatable    :: Smol_B_l(:)
  END TYPE Vec_smol_t

  TYPE :: psi_t
    TYPE (Basis_t),    pointer     :: Basis => Null()
    real (kind=Rk),    allocatable :: RVec(:)
    complex (kind=Rk), allocatable :: CVec(:)
    Real(kind=Rk), allocatable     :: Smol_G(:)
    Real(kind=Rk), allocatable     :: Smol_B(:)
    TYPE (Vec_smol_t), allocatable :: Vec_smol(:)
  END TYPE psi_t

  public :: psi_t,write_psi,init_psi,dealloc_psi
  ! operation on psi has to be defined: psi=psi1, psi1+psi2, psi=psi1*cte ...
contains
  SUBROUTINE init_psi(psi,Basis,cplx)
  USE Basis_m
    TYPE (psi_t),    intent(inout)      :: psi
    TYPE (Basis_t), intent(in), target :: Basis
    logical,        intent(in)         :: cplx

    CALL dealloc_psi(psi)


    IF (.NOT. Basis_IS_allocated(Basis)) STOP 'ERROR in Set_Op: the Basis is not initialized'


    IF (Basis%nb < 1) STOP 'ERROR in init_psi: Basis%nb < 1!'

    psi%Basis => Basis

    IF (cplx) THEN
      allocate(psi%CVec(Basis%nb))
    ELSE
      allocate(psi%RVec(Basis%nb))
    END IF

  END SUBROUTINE init_psi

  SUBROUTINE dealloc_psi(psi)
    TYPE(psi_t), intent(inout) :: psi

    nullify(psi%Basis)

    IF (allocated(psi%RVec)) THEN
      deallocate(psi%RVec)
    END IF

    IF (allocated(psi%CVec)) THEN
      deallocate(psi%CVec)
    END IF

  END SUBROUTINE dealloc_psi

  SUBROUTINE write_psi(psi)
    TYPE(psi_t), intent(in) :: psi

    IF (associated(psi%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to psi.'
    END IF

    IF (allocated(psi%RVec)) THEN
      write(out_unitp,*) 'Writing psi (real):'
      write(out_unitp,*) psi%RVec
      write(out_unitp,*) 'END Writing psi'
    END IF
    IF (allocated(psi%CVec)) THEN
      write(out_unitp,*) 'Writing psi (complex):'
      write(out_unitp,*) psi%CVec
      write(out_unitp,*) 'END Writting psi'
    END IF
  END SUBROUTINE write_psi
end module psi_m
