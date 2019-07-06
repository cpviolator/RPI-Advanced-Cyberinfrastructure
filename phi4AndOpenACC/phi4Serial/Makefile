#===========================================================
CXX=g++
CXXFLAGS = -std=c++11 -O3
#============================================================

phi4Serial: phi4Serial.o
	${CXX} ${CXXFLAGS} -o phi4Serial phi4Serial.o

phi4Serial.o: phi4Serial.cpp Makefile
	${CXX} ${CXXFLAGS} -c phi4Serial.cpp

clean:
	rm -f phi4Serial phi4Serial.o
