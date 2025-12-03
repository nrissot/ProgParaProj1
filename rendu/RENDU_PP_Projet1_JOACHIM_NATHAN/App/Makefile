# Compiler and flags
CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot

# Build any program from .cpp and run it immediately
.SUFFIXES:

arn_500_4 : main_arn run_arn_500_4 cleano

arn_500_16 : main_arn run_arn_500_16 cleano

arn_2k_4 : main_arn run_arn_2k_4 cleano

arn_2k_16 : main_arn run_arn_2k_16 cleano

# Pattern rule: build executable from .o
main_arn: main_arn.o ForARN.o FuncP.o Utils.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

run_arn_500_4:
	mpirun -np 4 --oversubscribe ./main_arn "exemples/dataset_500seq.fa" 4 500

run_arn_500_16:
	mpirun -np 16 --oversubscribe ./main_arn "exemples/dataset_500seq.fa" 4 500

run_arn_2k_4:
	mpirun -np 4 --oversubscribe ./main_arn "exemples/dataset_2000seq.fa" 4 2000

run_arn_2k_16:
	mpirun -np 16 --oversubscribe ./main_arn "exemples/dataset_2000seq.fa" 4 2000


# anti-frustration aliases
arn_2000_4 : arn_2k_4
arn_2000_16 : arn_2k_16
run_arn_2000_4: run_arn_2k_4
run_arn_2000_16: run_arn_2k_16

# Pattern rule: build object file from .cpp
%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Clean
cleano:
	rm -f *.o *~

clean:
	rm -f *.o *~ main_arn