CXX=mpicxx
CXXFLAGS= -std=c++11 -Wall -O2
LDLIBS=-lboost_program_options
TARGET=compile
default: $(TARGET)
all: $(TARGET)

Model.o: Model.cpp Model.h

Burgers.o: Burgers.cpp Burgers.h 

MPI_functions.o: MPI_functions.cpp MPI_functions.h 

burgers_main.o: burgers_main.cpp 

$(TARGET): burgers_main.cpp Burgers.o Model.o MPI_functions.o
	mpicxx burgers_main.cpp Burgers.o Model.o MPI_functions.o -o compile 

#Four serial cases, the parameters are: Px X Py, test case, number of points 
# in x (Nx), number of points in y (Ny)
diff: diff
	time mpiexec -np 1 ./$(TARGET) 1 100 100 100

advx: advx
	time mpiexec -np 1 ./$(TARGET) 2 100 100 100

advy: advy
	time mpiexec -np 1 ./$(TARGET) 3 100 100 100

burg: burg
	time mpiexec -np 1 ./$(TARGET) 4 100 100 100

#Four parallel cases, the parameters are: Px X Py, test case, number of points 
# in x (Nx), number of points in y (Ny)

diffp: diffp
	time mpiexec -np 2 ./$(TARGET) 1 100 100 100

advxp: advxp
	time mpiexec -np 2 ./$(TARGET) 2 100 100 100

advyp: advyp
	time mpiexec -np 2 ./$(TARGET) 3 100 100 100

burgp: burgp
	time mpiexec -np 2 ./$(TARGET) 4 100 100 100

#Target to clean the folder

.PHONY: clean 

clean:
	-rm -f *.o 

