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
PROGRAM TI_Schrod
  USE NumParameters_m
  USE Basis_m
  USE psi_m
  USE op_m
  USE NDindex_m
  USE Molec_m

  IMPLICIT NONE

  TYPE (Basis_t), target      :: Basis
  TYPE (psi_t)                :: psi,Hpsi
  TYPE (op_t)                 :: Op
!  TYPE (NDindex_t)            :: NDindex
  TYPE (Molec_t)              :: Molec
 !==============================================================================
 ! for QML
  REAL(kind=Rk), ALLOCATABLE  :: QQML(:)
  REAL(kind=Rk), ALLOCATABLE  :: Q(:)
  REAL(kind=Rk), ALLOCATABLE  :: Mat_V(:,:)
  REAL(kind=Rk)               :: V,Calc_pot1
  integer                     :: ndim,nsurf,option
  logical                     :: adiabatic
  character (len=16)          :: pot_name

  CALL Read_Molec(Molec,in_unitp)

  STOP
  !===================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  CALL Read_Basis(Basis,nio=in_unitp)

  CALL Set_op(Op,Basis,Molec) ! to be change
  CALL Make_Mat_OP(Op)
  CALL Diago_Op(Op)

  write(out_unitp,*) 'deallocation'

  CALL dealloc_Op(OP)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TI_Schrod
