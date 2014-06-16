#Makefile

HEADERS= geomap.h

OBJECTS= geomap.o

%.o : %.cpp $(HEADERS)
	gcc -ggdb -c -O3 -Wall $<

all : matchMesh

matchMesh: matchMesh.o $(OBJECTS)
	gcc -o matchMesh -Wall -L. matchMesh.o $(OBJECTS) -lm -lamrg -lstdc++

clean :
	rm -f *.o matchMesh
