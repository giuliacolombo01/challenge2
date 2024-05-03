CC = /u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++
CFLAGS = -I. -Ishared-folder/try -Wall -Werror

SRCS = main.cpp Matrix.hpp Matrix_more.hpp chrono.hpp
OBJS = $(SRCS:.cpp=.o)

all: main

main: $(OBJS)
	$(CXX) $(CFLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm *.o main
