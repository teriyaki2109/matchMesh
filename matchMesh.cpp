// *** PLEASE REFRAIN FROM DISTRIBUTING THIS BETA VERSION ***
// *** USE IT FOR RESEARCH PURPOSE ***
// Tony TUNG (tonytung.org)- 2014 (c)
// Surface Alignment by Geodesic Mapping [Tung et al., PAMI14] & [Tung et al., CVPR10]
// Use provided Makefile for compilation and don't forget to create an ./output/ folder before running the executable.
// Example of use: ./matchMesh toto/mesh001.off toto/mesh002.off 

#include <string>
#include <iostream>
#include "geomap.h"

int main(int argc,char **argv) 
{
	cout << "\n** Dynamic Surface Alignment using Geodesic Mapping **\n";
	if (argc < 3) return 0;

	//input files
	char *iFile=new char[100]; sprintf(iFile, "%s", argv[1]);
	char *jFile=new char[100]; sprintf(jFile, "%s", argv[2]); 
		
	match lMatch;
	lMatch.load(iFile, jFile);
	lMatch.computeExtremalPoints();
	lMatch.reorderReferences();
	lMatch.addReferenceIteratively();
	lMatch.matchMesh();
	
	return 0;
}




