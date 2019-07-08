#============================================================
TARGET  = 2D-Wilson
SOURCES = main.cpp
OBJS    = main.o
INC_PATH=../include

CXX=g++
CXXFLAGS = -O3 -g -std=c++11 -I${INC_PATH}
#============================================================

all: ${TARGET}

${TARGET}: main.o
	$(CXX) ${CXXFLAGS} -o ${TARGET} main.o

main.o: main.cpp Makefile
	${CXX} ${CXXFLAGS} -c main.cpp

#============================================================

clean:
	rm -f ${TARGET} main.o 
