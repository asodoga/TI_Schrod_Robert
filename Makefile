QMLDIR=/home/robert/robert2/Librairie/QuantumModelLib
AD_dnSVMDIR=/home/robert/robert2/Librairie/QuantumModelLib
ARPACKDIR=/home/robert/robert2/Librairie/ARPACK_DML
EVRTDIR=/home/robert/robert2/Librairie/ElVibRot-TnumTana/obj/obj_gfortran
PotDIR=/home/robert/robert2/potentiel/pot
PotRDIR=/home/robert/robert2/potentiel/pot/obj

ADDIR=$QMLDIR/Ext_Lib/dnSVMLib


FC=gfortran
#FFLAGS=-fopenmp
#FFLAGS=-O0 -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -fopenmp
#FFLAGS=-O0 -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -fopenmp
FFLAGS= -O0 -g -J obj -I$(PotRDIR) -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -fopenmp
#FFLAGS=-O5 -g -fbacktrace -funroll-loops -ftree-vectorize -falign-loops=16 -fopenmp
OBJ_DIR=obj
MOD_DIR=obj
SRC_DIR=src
MAIN_DIR=app



#QMLDIR=/Users/lauvergn/git/QuantumModelLib

#AD_dnSVMDIR=/Users/lauvergn/git/QuantumModelLib

#ARPACKDIR=/Users/lauvergn/trav/ARPACK





EXT_LIB=$(QMLDIR)/libpot.a $(PotDIR)/libpot.a $(AD_dnSVMDIR)/libAD_dnSVM.a $(ARPACKDIR)/libarpack_Linux.a




MAIN=TI_Schrod

LIBSRC=  UtilLib_m.f90 diago_m.f90 UtilMath_m.f90 NumParameters_m.f90  Tana_m.f90 \
 Molec_m.f90 NDindex_m.f90 Basis_m.f90 psi_m.f90 Op_m.f90
OBJ0=${LIBSRC:.f90=.o}

OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ $(OBJ_DIR)/files: $(OBJ))


$(MAIN).x: $(OBJ_DIR)/$(MAIN).o $(OBJ) $(EXT_LIB)
	$(FC) $(FFLAGS) -o $(MAIN).x  $(OBJ_DIR)/$(MAIN).o $(OBJ) $(EXT_LIB) -llapack -lblas

$(OBJ_DIR)/$(MAIN).o: $(MAIN_DIR)/$(MAIN).f90
	$(FC) $(FFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -o $@ -c $<


clean:
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

#QML, AD_dnSVM, ARPACK
$(QMLDIR)/libpot.a:
	test -e $(QMLDIR)/libpot.a || exit 1
$(AD_dnSVMDIR)/libAD_dnSVM.a:
	test -e $(AD_dnSVMDIR)/libAD_dnSVM.a || exit 1
$(ARPACKDIR)/libarpack_Linux.a:
	test -e $(ARPACKDIR)/libarpack_Linux.a || exit 1

#dependencies
$(OBJ_DIR)/$(MAIN).o: $(OBJ)

$(OBJ_DIR)/Basis_m.o: $(OBJ_DIR)/UtilMath_m.o $(OBJ_DIR)/NDindex_m.o $(OBJ_DIR)/UtilLib_m.o \
	 $(OBJ_DIR)/diago_m.o $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/Molec_m.o: $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/diago_m.o $(OBJ_DIR)/NumParameters_m.o

$(OBJ_DIR)/psi_m.o: $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/Op_m.o: $(OBJ_DIR)/Molec_m.o $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/NDindex_m.o: $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/diago_m.o:   $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/UtilLib_m.o: $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/UtilMath_m.o: $(OBJ_DIR)/NumParameters_m.o
