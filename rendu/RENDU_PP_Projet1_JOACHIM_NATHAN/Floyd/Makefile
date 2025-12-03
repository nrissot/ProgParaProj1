# Compiler and flags
CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot

# Build any program from .cpp and run it immediately
.SUFFIXES:

small: main run_small cleano

medium: main run_medium cleano

big: main run_big cleano

# Pattern rule: build executable from .o
main: main.o ForGraph.o FuncP.o Utils.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

run_small:
	mpirun -np 4 --oversubscribe ./main "exemples/PetitExemple.dot"

run_medium:
	mpirun -np 4 --oversubscribe ./main "exemples/MoyenExemple.dot"

run_big:
	mpirun -np 16 --oversubscribe ./main "exemples/GrosExemple.dot"


# Pattern rule: build object file from .cpp
%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Clean
cleano:
	rm -f *.o *~

clean:
	rm -f *.o *~ main