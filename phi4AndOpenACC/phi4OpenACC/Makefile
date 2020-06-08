#============================================================
#Uncomment the accelerator type with which you wish to compile

ACCEL_TYPE = PGI-tesla
#ACCEL_TYPE = PGI-multicore
#ACCEL_TYPE = PGI
#ACCEL_TYPE = GPP

TARGET=phi4_exe
CXX = pgc++
ifeq ($(ACCEL_TYPE),GPP)
        CXX=g++
        CXXFLAGS += -std=c++11 -O3 -DGPP_BUILD
else ifeq ($(ACCEL_TYPE),PGI)
        CXXFLAGS += -fast -std=c++11
        TARGET=phi4OpenACC_exe
else ifeq ($(ACCEL_TYPE),PGI-tesla)
        CXXFLAGS += -std=c++11 -Minfo=accel -fast -O3 -acc -ta=tesla:fastmath -ta:tesla:fma -ta=tesla:cc70 -DGPU -vvvv
        TARGET=phi4OpenACC-tesla_exe
else ifeq ($(ACCEL_TYPE),PGI-multicore)
        CXXFLAGS += -std=c++11 -Minfo=accel -fast -acc -ta=multicore
        TARGET=phi4OpenACC-multicore_exe
endif
#============================================================

v0: ${TARGET}_v0
v1: ${TARGET}_v1
v2: ${TARGET}_v2
v3: ${TARGET}_v3

all: v0 v1 v2 v3

${TARGET}_v0: phi4OpenACC-v0.o
	${CXX} ${CXXFLAGS} -o ${TARGET}_v0 phi4OpenACC-v0.o

phi4OpenACC-v0.o: phi4OpenACC-v0.cpp Makefile	
	${CXX} ${CXXFLAGS} -c phi4OpenACC-v0.cpp  


${TARGET}_v1: phi4OpenACC-v1.o
	${CXX} ${CXXFLAGS} -o ${TARGET}_v1 phi4OpenACC-v1.o

phi4OpenACC-v1.o: phi4OpenACC-v1.cpp Makefile	
	${CXX} ${CXXFLAGS} -c phi4OpenACC-v1.cpp  


${TARGET}_v2: phi4OpenACC-v2.o
	${CXX} ${CXXFLAGS} -o ${TARGET}_v2 phi4OpenACC-v2.o

phi4OpenACC-v2.o: phi4OpenACC-v2.cpp Makefile	
	${CXX} ${CXXFLAGS} -c phi4OpenACC-v2.cpp  


${TARGET}_v3: phi4OpenACC-v3.o
	${CXX} ${CXXFLAGS} -o ${TARGET}_v3 phi4OpenACC-v3.o

phi4OpenACC-v3.o: phi4OpenACC-v3.cpp Makefile	
	${CXX} ${CXXFLAGS} -c phi4OpenACC-v3.cpp  

clean:
	rm -f *_exe_v* phi4OpenACC*.o
