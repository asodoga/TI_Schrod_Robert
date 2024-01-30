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
!The main program which uses the subroutines of the other modules to solve
!the multidimensional Schrodinger equation using the direct product to evaluate
!the multiple integrations and then using the Smolyak Scheme
!===============================================================================
PROGRAM TI_Schrod
  USE NumParameters_m
  USE Basis_m
  USE psi_m
  USE op_m
  USE pot_m
  !USE NDindex_m
  USE Molec_m

  IMPLICIT NONE

  TYPE (Basis_t)              :: Basis
  TYPE (psi_t)                :: psi,Hpsi
  TYPE (op_t)                 :: Op
  TYPE (Molec_t)              :: Molec
  Real(kind=Rk), allocatable          :: Q(:),F2(:,:),F1(:)
  !!!!!!!!!!!!!!!!!test!!!!potlib!!!!!!!!!!!!!!!!!!!!!!!!!

  !====================================================================================
  !Allows you to read the information on the library used or not and the information
  !on the potential V of the systems being studied
  CALL Read_Molec(Molec,in_unitp)
!Stop 'Robert'
  !====================================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  !CALL Read_Basis_old(Basis,nio=in_unitp)

  !CALL Read_Basis_old1(Basis,nio=in_unitp)
  Call Read_Construct_Basis(Basis,nio=in_unitp)
  
  !!====================================================================================
  !The transfer of previously read information to the rest of the program
  CALL Set_op(Op,Basis,Molec) ! to be change
  !Stop 'Robert set_op'
  !Test smolyak!!!!!


  !====================================================================================
  !The action of the Hamiltonian operator on the wave packet and
  !the construction of the Hamiltonian matrix
  !CALL Make_Mat_OP2(Op)
  !CALL Make_Mat_OP_gen(Op,psi)
  CALL Make_mat(Op)
!
  !====================================================================================
  !The diagonalization of the Hamiltonian matrix
  CALL Diago_Op(Op)
  !STOP 'RObert A'
  !====================================================================================
  !Deallocating arrays allocated during the program
  write(out_unitp,*) 'deallocation'
  CALL dealloc_Op(OP)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TI_Schrod
