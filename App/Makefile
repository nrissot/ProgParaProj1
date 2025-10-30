# Compiler and flags
CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot

# Build any program from .cpp and run it immediately
.SUFFIXES:

all: main run cleano

# Pattern rule: build executable from .o
main: main.o ForGraph.o FuncP.o Utils.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

run:
	mpirun -np 16 --oversubscribe ./main "exemples/test.dot" 6

# Pattern rule: build object file from .cpp
%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Clean
cleano:
	rm -f *.o *~

clean:
	rm -f *.o *~ main