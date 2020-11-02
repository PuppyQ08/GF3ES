EXENAME = scf
OBJS = main.o scf.o

CXX = g++
CXXFLAGS = -std=c++1y  -c -g -O2 -Wall -Wextra -pedantic
LD = g++
LDFLAGS = -std=c++1y -lpthread -lm  -larmadillo

all : $(EXENAME)

$(EXENAME) : $(OBJS)
	$(LD) $(OBJS) $(LDFLAGS) -o $(EXENAME)

main.o : main.cpp scf.h
	$(CXX) $(CXXFLAGS) main.cpp

scf.o : scf.cpp scf.h
	$(CXX) $(CXXFLAGS) scf.cpp

clean :
	-rm -f *.o $(EXENAME)
