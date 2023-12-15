CXX = g++
CXXFLAGS = -Wall -Wextra -pthread -std=c++17

all: main
################################################
main: main.cc
	$(CXX) $(CXXFLAGS) -o main main.cc

#####################################################
.Phony: all clean
clean:
	rm -f *.o main