#ifndef _IMG_LOADER_H_
#define _IMG_LOADER_H_

#include "framework.h"
#include <memory>
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

std::unique_ptr<image<rgb> > LoadGDIPlusToRGB(const char* filename, imgFormat* format = 0);
std::unique_ptr<image<rgb> > LoadImageGenericRGB(const char* filename, imgFormat* format = 0);
std::unique_ptr <image<double> > LoadImageGenericMono(const char* filename,imgFormat* format=0);

void SaveImageGeneric(const image<rgb>* img,const char* filename,imgFormat format, int jpegQuality = 75);
void SaveImageGeneric(const image<double>* img,const char* filename,imgFormat format);

}

#endif //_IMG_LOADER_H_
