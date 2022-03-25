PROGRAM TI_Schrod
  USE NumParameters_m
  USE Basis_m
  USE psi_m
  USE op_m
  IMPLICIT NONE

  TYPE (Basis_t), target :: Basis
  TYPE (psi_t)           :: psi,Hpsi
  TYPE (op_t)            :: Op

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
!  CALL TEST_OpPsi_grid(Basis)

  CALL Set_op(Op,Basis) ! to be change
  CALL Make_Mat_OP(Op)
  CALL Diago_Op(Op)
!Stop '34'
  !CALL calc_OpPsi(H,psi,Hpsi)
!  CALL Write_psi(Hpsi)


  write(out_unitp,*) 'deallocation'
  CALL dealloc_Op(OP)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TI_Schrod
