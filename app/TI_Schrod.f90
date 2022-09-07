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
  !integer                     :: ndim,nsurf,option
  !logical                     :: adiabatic
  !character (len=16)          :: pot_name

!  ndim      = 0
!  nsurf     = 0
!  pot_name  = 'read_model'
!  adiabatic = .FALSE.
  !option    = 0
!  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
  CALL Read_Molec(Molec,in_unitp)

  ALLocate(Mat_V(1,1))
  ALLocate(QQML(3))
  QQML(:)= [1.6525859462_Rk,2.4849898659_RK,ZERO]


  CALL sub_Qmodel_V(Mat_V,QQML)
  Write(*,*) 'Mat_V',Mat_V
  ALLocate(Q(3))
  Q(:)= [2.4849898659_RK,2.4849898659_RK,1.6525859462_Rk]
  V = Calc_pot(Q)
  Write(*,*) 'V',V
  CALL Calc_potsub(Calc_pot1,Q,Molec)
  !STOP
  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  CALL Read_Basis(Basis,nio=in_unitp)
  !====================================================================

  !write(out_unitp,*) 'Initialization of a real psi'

  !CALL init_psi(psi,Basis,cplx=.FALSE.) ! to be changed
  !psi%RVec(:) = ONE
  !CALL Write_psi(psi)

!  write(out_unitp,*) 'Initialization of a complex psi'
!  CALL init_psi(psi,Basis,cplx=.TRUE.) ! to be changed
!  psi%CVec(:) = CONE
!  CALL Write_psi(psi)
!  write(out_unitp,*) ' | H | Psi > calculation'
  !Test Robert
  !CALL Test_Passage(Basis)
  !Call Calc_dngg_grid(Basis)


  CALL Set_op(Op,Basis) ! to be change
  CALL Make_Mat_OP(Op,Molec)
  CALL Diago_Op(Op,Molec)
!  CALL TEST_OpPsi_grid(Op)
!Stop '34'
  !CALL calc_OpPsi(H,psi,Hpsi)
!  CALL Write_psi(Hpsi)


  write(out_unitp,*) 'deallocation'

  CALL dealloc_Op(OP)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TI_Schrod
