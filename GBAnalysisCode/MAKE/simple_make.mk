srcDIR = $(DIR)

srcOBJGQS = \
Computers/cijkl.o \
Resources/Resources.o \
Resources/Exception.o \
Resources/index_map.o \
Resources/RNG_taus.o

#Compiler
CPP=icpc
CC=icc
FF=ifort

#Compiler flags
I_ACCURACY_FLAGS = -fp-model precise -prec-div -prec-sqrt
I_OPTIMIZATION_FLAGS = -O3 -xhost -ipo
FLAGS = $(I_OPTIMIZATION_FLAGS) $(I_ACCURACY_FLAGS) #-opt-report-file opt_report.txt
FFLAGS = $(FLAGS)
CFLAGS = $(FLAGS) -wr1125 #-wr21 -wr279 -wr1125 #-wr418 
LinkFLAGS = $(FLAGS)




INCLUDE = \
-I$(srcDIR) \
-I/data1/jcode/local/include 
#-I/data1/jamming/cpp/arpack++/include \
-I/data1/jcode/local/include/SuiteSparse \
-I/data1/jcode/local/include/netcdf4 \

LIBRARY = \
-L/data1/jamming/lib \
-L/data1/jcode/local/lib \
-L/data1/jcode/local/lib/static \
-L/usr/global/intel/Compiler/11.1/064/lib/intel64

SuiteSparseLINK = #-lumfpack -lamd -lcholmod -lcolamd -lblas
netCDFLINK = #-lnetcdf_c++4 -lnetcdf
hdf5LINK = #-lhdf5_hl -lhdf5 -lz
intelLINK = -lifcore -limf -lm
arpackLINK = #-larpack
LINK = $(SuiteSparseLINK)  -lgfortran $(arpackLINK) $(netCDFLINK) $(hdf5LINK) $(intelLINK) $(SuiteSparseLINK)


FRULE = $(FF)  $(FFLAGS) $(INCLUDE) -c -o $@ $<
CRULE = $(CPP) $(CFLAGS) $(INCLUDE) -c -o $@ $<
ORULE = $(CPP) $(CFLAGS) -o $@ $(OBJGQS) $(LIBRARY) $(LINK)


StandardDependencies = \
	$(srcDIR)/Boundaries/*.h \
	$(srcDIR)/Computers/*.h \
	$(srcDIR)/Database/*.h \
	$(srcDIR)/Minimization/*.h \
	$(srcDIR)/Potentials/*.h \
	$(srcDIR)/Resources/*.h \
	$(srcDIR)/State/*.h


#    $@  means "the target"
#    $<  means "whatever the dependencies are"



#####g++ -I/data1/jamming/cpp/arpack++/include -I/usr/global/netcdf-4.1.1-i11/include Resources/Resources.cpp Resources/Exception.cpp Test.cpp -L/data1/jamming/lib -L/data1/jamming/lib/lib_suitesparse -lamd -lcholmod -lcolamd -lccolamd -lcamd -lumfpack -lgfortran -larpack -L/usr/global/netcdf-4.1.1-i11/lib -lnetcdf_c++ -lnetcdf  -L/usr/global/hdf5-1.8.5-patch1-i11/lib -lhdf5_hl -lhdf5 -lz -lm -O3 -msse -msse2 -msse3 -o build/test.out
