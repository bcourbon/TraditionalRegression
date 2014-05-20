#makefile


CC = g++

#UCFLAGS = -O0 -g3 -Wall -gstabs+
UCFLAGS = -O3 -fopenmp -Wall -gstabs+
#UCFLAGS = -O0 -g3 -fopenmp -Wall -gstabs+

#RUCFLAGS = -pthread -m64 -I/afs/cern.ch/sw/lcg/external/root/head/slc4_amd64_gcc34/root/include
#LIBS = -L/afs/cern.ch/sw/lcg/external/root/head/slc4_amd64_gcc34/root/lib -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -pthread -lm -ldl -rdynamic
#GLIBS = -L/afs/cern.ch/sw/lcg/external/root/head/slc4_amd64_gcc34/root/lib -lCore -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGui -pthread -lm -ldl -rdynamic

RUCFLAGS := $(shell root-config --cflags) -I./include/
LIBS := -lgomp $(shell root-config --libs) -lTreePlayer -lTMVA -lRooFit -lRooFitCore 
GLIBS := $(shell root-config --glibs)

vpath %.cpp ./src

SRCPP = main.cpp\
Utilities.cpp\
GBRApply.cpp\
GBREvent.cpp\
GBRForest.cpp\
GBRTrainer.cpp\
GBRTree.cpp\
TMVAMaker.cpp\
RegressionTest.cpp\
GBRMaker.cpp\
ParReader.cpp\
RegressionManager.cpp\
SmearingCorrection.cpp\
ErrorCorrection.cpp\
TrackMomentumCorrection.cpp\
VariableCorrectionApply.cpp

#OBJCPP = $(SRCPP:.cpp=.o)
OBJCPP = $(patsubst %.cpp,obj/%.o,$(SRCPP))


all : regression.exe obj/libDictionary_C.so

obj/%.o : %.cpp
	@echo "> compiling $*"
	@mkdir -p obj/
	@$(CC) -c $< $(UCFLAGS) $(RUCFLAGS) -o $@

regression.exe : $(OBJCPP)
	@echo "> linking"
	@$(CC) $^ $(ACLIBS) $(LIBS) $(GLIBS) -o $@

clean:
	@echo "> Cleaning object files"
	@rm -f obj/*.o
        
cleanall: clean
	@echo "> Cleaning dictionary"
	@rm -f obj/libDictionary_C.so
	@echo "> Cleaning executable"
	@rm -f regression.exe

obj/libDictionary_C.so: ./include/libDictionary.C
	@echo "> Generating dictionary"
	@cd include && root -b -q libDictionary.C++
	@mv ./include/libDictionary_C.so ./obj/






