#ifndef _IMG_LOADER_H_
#define _IMG_LOADER_H_

#include "framework.h"

#include "image.h"

namespace img
{

typedef enum
{
	INVALID=0,
	PNM,
	BMP,
	GIF,
	JPEG, 
	PNG, 
	TIFF
}imgFormat;

void Startup();
void Shutdown();

image<rgb>* LoadGDIPlusToRGB(const char* filename, imgFormat* format = 0);
image<rgb>* LoadImageGenericRGB(const char* filename, imgFormat* format = 0);
image<double>* LoadImageGenericMono(const char* filename,imgFormat* format=0);

void SaveImageGeneric(const image<rgb>* img,const char* filename,imgFormat format, int jpegQuality = 75);
void SaveImageGeneric(const image<double>* img,const char* filename,imgFormat format);

}

#endif //_IMG_LOADER_H_
