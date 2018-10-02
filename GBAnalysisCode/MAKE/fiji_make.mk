

INCLUDE = \
-I$(srcDIR) \
-I/data1/jamming/cpp/arpack++/include \
-I/data1/jcode/local/include/SuiteSparse \
-I/usr/global/netcdf-4.1.1-i11/include \
-I/data1/jcode/local/include
#-I/usr/global/netcdf-4.1.1-i11/include \
-I/data0/home/cpgoodri/jmodes/cpp/netcdf-4.1.3/include \
-I/data1/jcode/local/include/netcdf4 \

LIBRARY = \
-L/data1/jamming/lib \
-L/data1/jcode/local/lib \
-L/usr/global/netcdf-4.1.1-i11/lib -lnetcdf_c++ \
-L/usr/global/hdf5-1.8.5-patch1-i11/lib \
-L/usr/global/intel/Compiler/11.1/064/lib/intel64
#-L/data1/jcode/local/lib/static \

SuiteSparseLINK = -lumfpack -lamd -lcholmod -lcolamd -lblas
netCDFLINK = -lnetcdf_c++ -lnetcdf
#netCDFLINK = -lnetcdf_c++4 -lnetcdf
hdf5LINK = -lhdf5_hl -lhdf5 -lz
intelLINK = -lifcore -limf -lm
LINK = $(SuiteSparseLINK)  -lgfortran -larpack $(netCDFLINK) $(hdf5LINK) $(intelLINK) $(SuiteSparseLINK)

