CC = g++
CXXFLAGS = -std=c++20
CPPFLAGS = -Wall -03 -I include

DOXYFILE = Doxyfile

SRCS = main.cpp Matrix.hpp Matrix_more.hpp chrono.hpp
OBJS = $(SRCS:.cpp=.o)

all: main

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm *.o main

doc:
	doxygen $(DOXYFILE)