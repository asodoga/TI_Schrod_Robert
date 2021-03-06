PROGRAM TI_Schrod
  USE NumParameters_m
  USE Basis_m
  USE psi_m
  USE op_m
  USE NDindex_m

  IMPLICIT NONE

  TYPE (Basis_t), target :: Basis
  TYPE (psi_t)           :: psi,Hpsi
  TYPE (op_t)            :: Op
!  TYPE (NDindex_t)       :: NDindex
 !==============================================================================
 ! for QML
  integer :: ndim,nsurf,option
  logical :: adiabatic
  character (len=16)                  :: pot_name

  ndim      = 0
  nsurf     = 0
  pot_name  = 'read_model'
  adiabatic = .FALSE.
  option    = 0
  CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
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
  CALL Make_Mat_OP(Op)
  CALL Diago_Op(Op)
!  CALL TEST_OpPsi_grid(Op)
!Stop '34'
  !CALL calc_OpPsi(H,psi,Hpsi)
!  CALL Write_psi(Hpsi)


  write(out_unitp,*) 'deallocation'

  CALL dealloc_Op(OP)
  CALL dealloc_psi(psi)
  CALL dealloc_psi(Hpsi)

END PROGRAM TI_Schrod
