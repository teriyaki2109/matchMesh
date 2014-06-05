#ifndef __BYTEORDER_H_
#define __BYTEORDER_H_

#ifdef _WIN32
#define LITTLE_ENDIAN 1234
#endif

#ifndef BYTE_ORDER
#if defined _LITTLE_ENDIAN || defined LITTLE_ENDIAN
#ifndef LITTLE_ENDIAN
#define LITTLE_ENDIAN 1234
#endif
#define BYTE_ORDER LITTLE_ENDIAN
#elif defined _BIG_ENDIAN || defined BIG_ENDIAN
#ifndef BIG_ENDIAN
#define BIG_ENDIAN 4321
#endif
#define BYTE_ORDER BIG_ENDIAN
#endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

size_t readOtherE(void * ptr,unsigned int sz,unsigned int cpt,FILE* f);

#endif
