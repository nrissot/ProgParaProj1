CXX = mpic++
CXXFLAGS = -std=c++11 -Wall -Wextra -pedantic -lcgraph -lgvc -lxdot -fopenmp

main: main.o needleman_wunsch.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

run:
	./main

clean:
	rm -f *.o *~ main