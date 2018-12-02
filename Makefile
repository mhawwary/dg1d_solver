
# Specifying some environmental variables for Linux, note this can be done in the shell prompt

COMP	= GCC
TECIO	= NO
CODE	= RELEASE
OPENMP	= NO

ifeq ($(CODE),RELEASE)
     bin_name = DG1DFlow.exe
endif

ifeq ($(CODE),DEBUG)
    bin_name = DG1DFlow_debug.exe
endif

# Specifing Standard Variables:
CXX	= g++ -std=gnu++11 #-pedantic-errors # c++ gcc compiler
CXXFLAGS=       # C++ compiler flags
LDLFLAGS=	# linker flags
CPPFLAGS=	# c/c++ preprocessor flags

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
DEBUG_dir = debug/
INC	= include/

vpath %.cpp src include
vpath %.c src
vpath %.o   obj
vpath %.h include src
vpath %.hpp include src

# Objects
OBJS	= $(OBJ)DG1DFlow.o $(OBJ)SimData.o $(OBJ)ExplicitTimeSolver.o $(OBJ)DGSolverAdvec.o $(OBJ)DGSolverDiffus.o $(OBJ)DGSolverAdvecDiffus.o $(OBJ)Faddeeva.o # objects 
INCLS	= 

# Compile

.PHONY: default help clean


default: all
help:	
	@echo 'help'

all: $(bin_name)

ifeq ($(CODE),DEBUG)
    $(bin_name): $(OBJS) 
	$(CXX) $(OPTS) -o $(DEBUG_dir)$@ $+
endif
ifeq ($(CODE),RELEASE)
    $(bin_name): $(OBJS) 
	$(CXX) $(OPTS) -o $(BIN)$@ $+
endif

$(OBJ)%.o : %.cpp 
	$(CXX) $(OPTS) -c -o $@ $<


$(OBJ)DG1DFlow.o:   DG1DFlow.cpp 
$(OBJ)SimData.o:   SimData.hpp SimData.cpp
$(OBJ)ExplicitTimeSolver.o: ExplicitTimeSolver.hpp ExplicitTimeSolver.cpp
$(OBJ)DGSolverAdvec.o: DGSolverAdvec.hpp
$(OBJ)DGSolverDiffus.o: DGSolverDiffus.hpp
$(OBJ)DGSolverAdvecDiffus.o: DGSolverAdvecDiffus.hpp
$(OBJ)Faddeeva.o: Faddeeva.hpp Faddeeva.cpp 

clean:
	rm -f ./$(OBJ)*.o ./$(BIN)*.exe  
	@echo  removing all object and executable files

clean_debug:
	rm -f ./$(OBJ)*.o ./$(DEBUG_dir)*.exe  
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


