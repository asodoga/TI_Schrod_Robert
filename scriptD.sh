#!/bin/bash
#QMLDIR=/home/robert/robert2/CODE_CORECT/QuantumModelLib
QMLDIR=/Users/lauvergn/git/QuantumModelLib
FPMDIR=build/gfortran*/TI_Schrod


ADDIR=$QMLDIR/Ext_Lib/dnSVMLib

Tana_file=TanaF2F1Vep_ClD2p.f90

#make Tana_m.f90
echo "module Tana_m
  implicit none
contains" > src/Tana_m.f90
cat Tana_files/$Tana_file >> src/Tana_m.f90
echo "END module Tana_m" >> src/Tana_m.f90
#end make Tana_m.f90


rm -r build
fpm build --flag  "-fopenmp"
gfortran -fopenmp -o TI_Schrod.x  $FPMDIR/app_TI_Schrod.f90.o $FPMDIR/libTI_Schrod.a $QMLDIR/libpot.a $ADDIR/libAD_dnSVM.a -lblas -llapack

./TI_Schrod.x < DAT_files/Dat_test_QML > Reread
