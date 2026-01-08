CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp

all: main cleano

2000: main run_2000 cleano

500: main run_500 cleano

20: main run_20 cleano

run_2000:
	mpirun -np 16 --oversubscribe ./main "arn/dataset_2000seq.fa" 4 2000

run_500:
	mpirun -np 4 --oversubscribe ./main "arn/dataset_500seq.fa" 4 500

run_20:
	mpirun -np 4 --oversubscribe ./main "arn/dataset_20seq.fa" 4 20

main: main.o needleman_wunsch.o ForARN.o FuncP.o Utils.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f *.o *~ main

cleano:
	rm -f *.o *~