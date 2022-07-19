module Molec_m
  USE NumParameters_m
  implicit none
  private

  !real(kind=Rk) :: mass = HALF
  real(kind=Rk) :: mass3(3) != HALF
  real(kind=Rk) :: massr1
  real(kind=Rk) :: massr2
  real(kind=Rk) :: mass = 1744.44536_RK
  !real(kind=Rk) :: mass = ONE
  real(kind=Rk) :: De = 0.2250_RK
  real(kind=Rk) :: Re = 1.73290_RK
  real(kind=Rk) :: alpha1 = 1.1741_RK

  public :: Calc_pot,mass,mass3,massr1,massr2

  !massr1= (mass3(1)*mass3(3))/(mass3(1)+mass3(3))
  !massr2= (mass3(1)*mass3(3))/(mass3(1)+mass3(3))

contains
  FUNCTION Calc_pot(Q)
    real(kind=Rk)             :: Calc_pot

    real(kind=Rk), intent(in) :: Q(:)

  !Calc_pot = HALF * dot_product( Q,Q)! 0.5*x^2
  !Calc_pot =  dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! x^2+x^4
  !Calc_pot =  -TEN * dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! -10x^2+x^4
  !Calc_pot = HALF * dot_product( Q,Q)
    !Calc_pot = De*dot_product((1-exp(-alpha1*(Q-Re))),(1-exp(-alpha1*(Q-Re))))
    Calc_pot = De*(1-exp(-alpha1*(Q(1)-Re)))**2
  END FUNCTION Calc_pot

end module Molec_m
