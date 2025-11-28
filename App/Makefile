# Compiler and flags
CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot

# Build any program from .cpp and run it immediately
.SUFFIXES:

all: main run cleano

arn: main_arn run_arn_2000 cleano

arn_500: main_arn run_arn_500 cleano

arn_20: main_arn run_arn_20 cleano

# Pattern rule: build executable from .o
main: main.o ForGraph.o FuncP.o Utils.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

main_arn: main_arn.o ForARN.o FuncP.o Utils.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

run:
	mpirun -np 16 --oversubscribe ./main "exemples/exported_graph.dot" 4

run_arn_2000:
	mpirun -np 16 --oversubscribe ./main_arn "exemples/arn/dataset_2000seq.fa" 4 2000

run_arn_500:
	mpirun -np 16 --oversubscribe ./main_arn "exemples/arn/dataset_500seq.fa" 4 500

run_arn_20:
	mpirun -np 16 --oversubscribe ./main_arn "exemples/arn/dataset_20seq.fa" 4 20

mem:
	valgrind mpirun -np 16 --oversubscribe ./main "exemples/exported_graph.dot" 4

# Pattern rule: build object file from .cpp
%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Clean
cleano:
	rm -f *.o *~

clean:
	rm -f *.o *~ main main_arn