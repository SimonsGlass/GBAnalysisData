DIR = /data0/home/cpgoodri/jcode/test_project/jsrc

name=DiagonalizationTutorial
obj=$(name).o



prjDIR = .
prjOBJGQS = \
$(obj)

include $(DIR)/std_make.mk

OBJGQS = \
$(patsubst %,$(srcDIR)/%,$(srcOBJGQS)) \
$(patsubst %,$(prjDIR)/%,$(prjOBJGQS))




.f.o: 
	$(FRULE)

.cpp.o: 
	$(CRULE)

$(name).out: $(OBJGQS)
	$(ORULE)

#if any header file is changed, the project file gets recompiled.
$(obj): $(StandardDependencies)





