CC = /u/sw/toolchains/gcc-glibc/11.2.0/prefix/bin/g++
CFLAGS = -I. -Ishared-folder/try -Wall -Werror

SRCS = main.cpp Matrix.h Matrix_more.h chrono.h
OBJS = $(SRCS:.cpp=.o)

all: main

main: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) main