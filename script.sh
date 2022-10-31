#!/bin/bash
QMLDIR=/home/robert/robert2/Librairie/QuantumModelLib
ARPACKDIR=/home/robert/robert2/Librairie/ARPACK_DML
EVRTDIR=/home/robert/robert2/Librairie/ElVibRot-TnumTana/obj/obj_gfortran
FPMDIR=build/gfortran*/TI_Schrod


ADDIR=$QMLDIR/Ext_Lib/dnSVMLib

#Tana_file=TanaF2F1Vep.f90
#Tana_file=TanaF2F1Vep_ClD2p.f90
#Tana_file=TanaF2F1Vep_ClHDp.f90
Tana_file=TanaF2F1Vep_37ClH2p.f90
#make Tana_m.f90
echo "module Tana_m
  implicit none
contains" > src/Tana_m.f90
cat Tana_files/$Tana_file >> src/Tana_m.f90
echo "END module Tana_m" >> src/Tana_m.f90
#end make Tana_m.f90
rm -r build
fpm build --flag  "-fopenmp"
echo "fpm done"
#gfortran -fopenmp -o TI_Schrod.x  $FPMDIR/app_TI_Schrod.f90.o $FPMDIR/libTI_Schrod.a $QMLDIR/libpot.a $ADDIR/libAD_dnSVM.a -lblas -llapack
gfortran -fopenmp -o TI_Schrod.x  $FPMDIR/app_TI_Schrod.f90.o $FPMDIR/libTI_Schrod.a $QMLDIR/libpot.a $ADDIR/libAD_dnSVM.a $ARPACKDIR/libarpack_Linux.a $EVRTDIR/libTnum.a -lblas -llapack
./TI_Schrod.x < DAT_files/Dat_test_QML_op6_37clh2p > arpacka
#./TI_Schrod.x < DAT_files/Dat_test_QML_ClH2+_Botschwina2 > 10nQML_ClH2+_Botschwina2

#./TI_Schrod.x < Resultats_ClH2p/Dat_test_QML_ClH2+_Botschwina > ClD2+_Botschwina1

#./TI_Schrod.x < Resultats_ClDDp/Dat_test_QML_op6 > restClHDp_QML_op6
#./TI_Schrod.x < Resultats_ClDDp/Dat_test_QML_ClH2+_Botschwina > ClD2+_Botschwina1

#./TI_Schrod.x < Resultats_ClHDp/Dat_test_QML_op6 > restClHDp_QML_op6
#./TI_Schrod.x < Resultats_ClHDp/Dat_test_QML_ClH2+_Botschwina > ClD2+_Botschwina1
