module Molec_m
  USE NumParameters_m
  implicit none
  private

  real(kind=Rk) :: mass = HALF
  !real(kind=Rk) :: mass = ONE

  public :: Calc_pot,mass

contains
  FUNCTION Calc_pot(Q)
    real(kind=Rk)             :: Calc_pot

    real(kind=Rk), intent(in) :: Q(:)

  !Calc_pot = HALF * dot_product( Q,Q)! 0.5*x^2
  !Calc_pot =  dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! x^2+x^4
    Calc_pot =  -TEN * dot_product( Q,Q) + dot_product( Q,Q) *dot_product( Q,Q)! -10x^2+x^4

  !Calc_pot = HALF * dot_product( Q,Q)

  END FUNCTION Calc_pot

end module Molec_m
