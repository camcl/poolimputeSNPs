## prophaser: Compilation with GDB for debugging
CC = g++
CFLAGS = -std=c++14 -Ofast -g3 -Wall -c -fopenmp -msse2 -mavx -DNDEBUG
LFLAGS = -fopenmp -g -o

SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=phase
LIBSTATGEN=/opt/conda/pkgs/libstatgen-*
EIGEN=/opt/conda/pkgs/eigen-*/include/eigen3

all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CC) $(LFLAGS) $@ $^  -L$(LIBSTATGEN) -lStatGen -lz
#	$(CC) $(LFLAGS) $@ $^  -L$(LIBSTATGEN) -lStatGen_debug -lz
%.o: %.cpp %.h
	$(CC) $(CFLAGS) -g $< -I $(LIBSTATGEN)/include/ -I $(EIGEN)


%.o: %.cpp
	 $(CC) $(CFLAGS) $< -I $(LIBSTATGEN)/include/  -I $(EIGEN)


clean:
	rm -f *.o

