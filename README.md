TI_Schrod_Robert is code in development to solve the schrodinger equation for large systems.
date: 02/02/2022
Copyright 2022 Robert AFANSOUNOUDJI [1]

     with contributions of:
     Komi SODOGA      [1]
     David Lauvergnat [2]

[1]:Laboratoire de Physique des Matériaux et des Composants à Semi-Conducteurs (LPMCS), UNIVERSITÉ DE LOME
[2]:Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay,

- General quantum dynamics code using curvilinear coordinates:
 1. Vibrational levels, intensities for soft molecular systems
 2. Optimization with the given set of curvilinear coordinates
 3. Calculation of the average value of positions of atoms of the molecule

- This code uses libraries like ElVibRot-TnumTana, QuantumModelLib, Lapack, arpack and Blas

- Originalities of this code:
   1. No built-in limitation in terms of number of degrees of freedom.
   2. the use of a numerical but exact kinetic energy operator with Tnum (ElVibRot-TnumTana) (Automatic Differentiation),    which allows great flexibility in the choice of curvilinear coordinates.
   2. the use of the Smolyak scheme, which avoids a basic set and grids of direct products under development..

   - List of constant names
    This list of names can be used with parameter. However, the following parameters can be used in this  list of names for information about the bases and potentials being investigated:

   * coord_type='CURV',
   * sym_type='SYM',
   * Model_type='QML', When using the potentials of the QML library,
   * V0 =  minimum value of the potential,
   * Q0 =  minimum value of the position  
   * pot_name = potential surface name
   * name ='Dp'or 'smoly' , Dp for direct product and smoly for smolyak
   * nb_basis = number of bases used   
   * name = 'herm 'or 'boxab' ,The name of the base used either herm for the base of harmonic   oscillators or boxab for the base of boxes (sine bases)  
   * nb = number of basic elements according to each degree of freedom
   * nq = number of grid points along each degree of freedom
   * SCALEQ =   scaling factor of the mathematical base used  



































Introduction

MODULE

 Basis_m :
  The derived type:
  Basis_t Which contains:
  * The integers like
  - nb_basis is number of bases used   
  - nb  is number of basic elements according to each degree of freedom        
  - nq  is number of grid points along each degree of freedom      
  * The character like
  - Basis_name is the name of the base  used
  * The reals like
  - x(:) are the positions of the grid points along each degree of freedom
  - w(:) are the weights of the grid points along each degree of freedom
  - d0gb(:,:)      are basis functions d0gb(nq,nb)
  - d1gb(:,:,:)    are basis functions d2gb(nq,nb,1)
  - d1gg(:,:,:)    are basis functions d2gg(nq,nq,1) , on the grid
  - d2gb(:,:,:,:)  are basis functions d2gb(nq,nb,1,1),on the grid
  - d2gg(:,:,:,:)  are basis functions d2gg(nq,nq,1,1),on the grid


 SUBROUTINES

 Read_Basis,

 Basis_IS_allocated,

 BasisTOGrid_Basis,

 GridTOBasis_Basis Test_Passage,

 Calc_dngg_grid,

 Basis_IS_allocatedtot,

 write_basis
