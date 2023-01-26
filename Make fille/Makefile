FC=gfortran
FFLAGS=-O3 -Wall -Wextra -fopenmp -J$(MOD_DIR)

OBJ_DIR=obj
MOD_DIR=obj
SRC_DIR=src
MAIN_DIR=app


#QMLDIR=/home/elprof/QuantumModelLib
QMLDIR=/home/robert/robert2/Librairie/QuantumModelLib
AD_dnSVMDIR=/home/robert/robert2/Librairie/ARPACK_DML
ARPACKDIR=/home/robert/robert2/Librairie/ElVibRot-TnumTana/obj/obj_gfortran

EXT_LIB=$(QMLDIR)/libpot.a $(AD_dnSVMDIR)/libAD_dnSVM.a $(ARPACKDIR)/libarpack_OSX.a

MAIN=TI_SCHROD

LIBSRC= NumParameters_m.f90 UtilLib_m.f90 diago_m.f90 \
 Molec_m.f90 NDindex_m.f90 Basis_m.f90 psi_m.f90 Op_m.f90 Tana_m.f90
OBJ0=${LIBSRC:.f90=.o}

OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ $(OBJ_DIR)/files: $(OBJ))


$(MAIN).x: $(OBJ_DIR)/$(MAIN).o $(OBJ) $(EXT_LIB)
	$(FC) $(FFLAGS) -o $(MAIN).x  $(OBJ_DIR)/$(MAIN).o $(OBJ) -lblas -llapack $(EXT_LIB)

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
$(ARPACKDIR)/libarpack_OSX.a:
	test -e $(ARPACKDIR)/libarpack_OSX.a || exit 1

#dependencies
$(OBJ_DIR)/$(MAIN).o: $(OBJ)

$(OBJ_DIR)/Basis_m.o: $(OBJ_DIR)/NDindex_m.o $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/diago_m.o $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/Molec_m.o: $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/diago_m.o $(OBJ_DIR)/NumParameters_m.o

$(OBJ_DIR)/psi_m.o: $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/Op_m.o: $(OBJ_DIR)/Molec_m.o $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/NDindex_m.o: $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/diago_m.o:   $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/UtilLib_m.o: $(OBJ_DIR)/NumParameters_m.o
