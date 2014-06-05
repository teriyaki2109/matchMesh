#ifndef __3DFILES_H_
#define __3DFILES_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "geombasic.h"
#include "object.h"
#include "list.h"
//#include "statusBar.h"
#include "byteOrder.h"
 
enum Format {VRML,OFF,TRI,STL,UNKNOWN};

void read3D(Object* &obj, char* name,bool* const &pointsToKeep, bool texture_check, bool acp_check);
void reloadObject(Object* &obj,char* name,bool* const &pointsToKeep, bool texture_check, bool acp_check);
void writeTRI(Object* const &obj,char* name,int resolution, bool texture_check);
void writeCOFF(Object* const &obj,char* name,int resolution);
Edge* addEdge(Object* obj,int indexPoint1,int indexPoint2,int* indexEdge,Face* const &_face);
#endif




