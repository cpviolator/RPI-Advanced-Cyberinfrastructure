#============================================================
TARGET = ads_heatbath
OBJS   = main.o util.o hyp_util.o graph.o
#============================================================
CXX = g++
CXXFLAGS = -O3 -g -std=c++11
#============================================================

all: ${TARGET}

main.o: main.cpp Makefile	
	${CXX} ${CXXFLAGS} -c main.cpp 

util.o: util.cpp 
	${CXX} ${CXXFLAGS} -c util.cpp

hyp_util.o: hyp_util.cpp 
	${CXX} ${CXXFLAGS} -c hyp_util.cpp

graph.o: graph.cpp 
	${CXX} ${CXXFLAGS} -c graph.cpp

${TARGET}: ${OBJS} 
	${CXX} $(CXXFLAGS) -o ${TARGET} ${OBJS}

#============================================================

clean:
	rm -f ${TARGET} ${OBJS} core*
