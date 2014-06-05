#ifndef __DATAITEM_H_
#define __DATAITEM_H_

#include <stdlib.h>
#include <string.h>

class DataItem{
 public:
  char title[64];
  float min,max;
  int length;
  float *data;
// DataItem(){};
  DataItem() :
	length(0),
	data(NULL)
	{
		memset(title,'\0',64);
	}
	
  ~DataItem(){ if (data) delete [] data;};
  DataItem(const char *tit,const float *dat,int len) :
    length(len)
    {
      strcpy(title,tit);
      //      length=len;
      data=new float[length];
      memcpy(data,dat,length*sizeof(float));
      min=(float)HUGE_VAL;
      max=(float)-HUGE_VAL;
      for(int i=0;i<length;i++){
        if(data[i]>max) max=data[i];
        if(data[i]<min) min=data[i];
      }
    };
  void write(FILE *fd)
    {
      fwrite(this,sizeof(DataItem)-sizeof(float*),1,fd);
      fwrite(data,sizeof(float)*length,1,fd);
    };
  void read(FILE *fd)
    {
      fread(this,sizeof(DataItem)-sizeof(float*),1,fd);
	  if (data) delete [] data;
      data=new float[length];
      fread(data,sizeof(float)*length,1,fd);
    };
};

#endif
