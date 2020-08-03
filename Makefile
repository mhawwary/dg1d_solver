
# Specifying some environmental variables for Linux, note this can be done in the shell prompt

COMP	= GCC
TECIO	= NO
CODE	= RELEASE
OPENMP	= NO
PLATFORM = LINUX

# Specifing Standard Variables:
ifeq ($(PLATFORM), LINUX)
	CXX = g++ #-pedantic-errors # c++ gcc compiler
	LDLFLAGS =  # linker flags
        BINNAME = dg1dflow 
else
	CXX = x86_64-w64-mingw32-g++ # c++ gcc compiler
	LDLFLAGS = -static-libgcc -static-libstdc++ # linker flags
        BINNAME = dg1dflow.exe
endif
CXXFLAGS= -std=gnu++11 # C++ compiler flags
CPPFLAGS=  # c/c++ preprocessor flags

OPTS	= 	# optimization flags and other options

# Includes

OPTS	+= -I include
#OPTS    += -I /home/mhawwary/Libraries/Eigen3.3.2/Eigen
#OPTS += -I /home/mhawwary/work/hpMusic/contrib/eigen/Eigen

ifeq ($(TECIO),YES)
	OPTS += -I $(TECIO_DIR)/tecsrc
endif

ifeq ($(CODE),RELEASE)
	ifeq ($(COMP),GCC)
		OPTS	+= -O3 
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif

ifeq ($(OPENMP),YES)	
	OPTS	+=  -lgomp -fopenmp 
endif

ifeq ($(CODE),DEBUG)
	ifeq ($(COMP),GCC)
		OPTS	+= -fbuiltin -g -Wall #-Werror
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif


# Source

SRC	= src/
OBJ	= obj/
BIN	= bin/
INC	= include/
$(shell mkdir obj)
$(shell mkdir bin)

vpath %.cpp src include
vpath %.c src
vpath %.o   obj
vpath %.h include src
vpath %.hpp include src

# Objects
OBJS	= $(OBJ)DG1DFlow.o $(OBJ)SimData.o $(OBJ)ExplicitTimeSolver.o $(OBJ)DGSolverAdvec.o $(OBJ)DGSolverDiffus.o $(OBJ)DGSolverAdvecDiffus.o $(OBJ)general_tools.o $(OBJ)Faddeeva.o # objects 
INCLS	= 

# Compile

.PHONY: default help clean


default: all
help:	
	@echo 'help'

all: $(BINNAME)
$(BINNAME): $(OBJS)
	$(CXX) $(OPTS) -o $(BIN)$@ $+

$(OBJ)%.o : %.cpp 
	$(CXX) $(CXXFLAGS) $(OPTS) -c -o $@ $<
$(OBJ)%.o : %.c 
	$(CXX) $(CXXFLAGS) $(OPTS) -c -o $@ $<


$(OBJ)DG1DFlow.o:   DG1DFlow.cpp 
$(OBJ)SimData.o:   SimData.hpp SimData.cpp
$(OBJ)ExplicitTimeSolver.o: ExplicitTimeSolver.hpp ExplicitTimeSolver.cpp
$(OBJ)DGSolverAdvec.o: DGSolverAdvec.hpp
$(OBJ)DGSolverDiffus.o: DGSolverDiffus.hpp
$(OBJ)DGSolverAdvecDiffus.o: DGSolverAdvecDiffus.hpp
$(OBJ)Faddeeva.o: Faddeeva.hpp Faddeeva.cpp 
$(OBJ)general_tools.o: general_tools.cpp

clean:
	rm -f ./$(OBJ)*.o ./$(BIN)*.exe ./$(BIN)$(BINNAME) 
	@echo  removing all object and executable files

clean_temp:
	rm -f ./input/*.*~ ./$(OBJ)*.*~ ./$(BIN)*.*~ ./$(SRC)*.*~ ./$(INC)*.*~ ./*.*~ ./*~
	rm -f ./Results/*.*~ ./Results/*/*.*~ ./Results/*/*/*.*~ ./Results/*/*/*/*.*~ 
	

plot:
	python python_tools/DGsolplot_reader.py -f ./input/python_input.in 

plot_compare:
	python python_tools/DGplot_compare_Beta.py -f ./input/python_input_compare_Beta.in 

plot_errors:
	python python_tools/DGplot_errors.py -f ./input/python_input_errors.in 

plot_error_analys:
	python python_tools/DGplot_error_analysis.py -f ./input/python_input_error_analysis.in 


