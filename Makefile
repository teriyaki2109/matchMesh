#Makefile

HEADERS= geomap.h

OBJECTS= geomap.o

%.o : %.cpp $(HEADERS)
	gcc -ggdb -c -O3 -Wall $<

all : matchMesh

matchMesh: matchMesh.o $(OBJECTS)
	gcc -o matchMesh matchMesh.o $(OBJECTS) -lm -Wall -lstdc++ -L. -lamrg

clean :
	rm -f *.o matchMesh
