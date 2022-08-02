#!/bin/bash
QMLDIR=/home/robert/robert2/CODE_CORECT/QuantumModelLib
FPMDIR=build/gfortran*/TI_Schrod


ADDIR=$QMLDIR/Ext_Lib/dnSVMLib

rm -r build
fpm build --flag  "-fopenmp"
gfortran -fopenmp -o TI_Schrod.x  $FPMDIR/app_TI_Schrod.f90.o $FPMDIR/libTI_Schrod.a $QMLDIR/libpot.a $ADDIR/libAD_dnSVM.a -lblas -llapack
#./TI_Schrod.x <DAT_files/dat_dia> resultat
./TI_Schrod.x < DAT_files/Dat_test_QML > Resultat
