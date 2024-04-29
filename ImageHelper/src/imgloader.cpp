#include <windows.h>
#include <gdiplus.h>
#include <stdio.h>
#include <algorithm>
#include "../inc/imgloader.h"
#include "../inc/pnmfile.h"
#include "../inc/misc.h"

using namespace pnm;

namespace img
{

static ULONG_PTR s_Token = 0UL;

void Startup()
{
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	Gdiplus::GdiplusStartup(&s_Token,&gdiplusStartupInput,NULL);
}

void Shutdown()
{
	Gdiplus::GdiplusShutdown(s_Token);
}
	
static int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
{
	UINT  num = 0;          // number of image encoders
	UINT  size = 0;         // size of the image encoder array in bytes

	Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;

	Gdiplus::GetImageEncodersSize(&num, &size);
	if(size == 0)return -1;  // Failure

	pImageCodecInfo = (Gdiplus::ImageCodecInfo*)(malloc(size));
	if(pImageCodecInfo == NULL)return -1;  // Failure

	Gdiplus::GetImageEncoders(num, size, pImageCodecInfo);

	for(UINT j = 0; j < num; ++j)
	{
#pragma warning( suppress: 6385 ) // static analysis doesn't understand sizes of the data
		if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 )
		{
			*pClsid = pImageCodecInfo[j].Clsid;
			free(pImageCodecInfo);
			return j;  // Success
		}    
	}

	free(pImageCodecInfo);
	return -1;  // Failure
}

static imgFormat GetFormat(Gdiplus::Image* gdiimage)
{
	GUID guid;
	gdiimage->GetRawFormat(&guid);
	if(guid==Gdiplus::ImageFormatBMP)
	{
		return BMP;
	}
	else if(guid==Gdiplus::ImageFormatGIF)
	{
		return GIF;
	}
	else if(guid==Gdiplus::ImageFormatJPEG)
	{
		return JPEG;
	}
	else if(guid==Gdiplus::ImageFormatPNG)
	{
		return PNG;
	}
	else if(guid==Gdiplus::ImageFormatTIFF)
	{
		return TIFF;
	}
	return INVALID;
}

image<rgb>* LoadGDIPlusToRGB(const char* filename, imgFormat* format) {
	ULONG_PTR gdiplusToken;
	if (s_Token == 0UL) {
		Gdiplus::GdiplusStartupInput gdiplusStartupInput;
		Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);
	}
	else {
		gdiplusToken = s_Token;
	}

	int nLen = MultiByteToWideChar(CP_ACP, 0, filename, -1, NULL, NULL);
	LPWSTR lpszW = new WCHAR[nLen];
	MultiByteToWideChar(CP_ACP, 0, filename, -1, lpszW, nLen);
	Gdiplus::Image* gdiimage = new Gdiplus::Image(lpszW);
	delete[] lpszW;
	if (gdiimage->GetLastStatus() != Gdiplus::Ok)
	{
		delete gdiimage;
		if (s_Token == 0UL) {
			Gdiplus::GdiplusShutdown(gdiplusToken);
		}
		return 0;
	}
	if (format)
	{
		*format = GetFormat(gdiimage);
	}
	int width = gdiimage->GetWidth();
	int height = gdiimage->GetHeight();
	Gdiplus::Bitmap* tempBMP = new Gdiplus::Bitmap(width, height, PixelFormat32bppPARGB);
	tempBMP->SetResolution(gdiimage->GetHorizontalResolution(), gdiimage->GetVerticalResolution());
	Gdiplus::Graphics* graphics = new Gdiplus::Graphics(tempBMP);
	graphics->DrawImage(gdiimage, 0, 0);
	delete graphics;
	delete gdiimage;
	Gdiplus::BitmapData bitmapData;
	Gdiplus::Rect fullRect(0, 0, width, height);
	tempBMP->LockBits(&fullRect, Gdiplus::ImageLockModeRead, PixelFormat32bppPARGB, &bitmapData);
	unsigned int* pRaw = (unsigned int*)bitmapData.Scan0;
	image<rgb>* output = new image<rgb>(width, height, false);
	for (int x = 0; x < width; x++)
	{
		for (int y = 0; y < height; y++)
		{
			unsigned int curColor = pRaw[y * bitmapData.Stride / 4 + x];
			uchar b = curColor & 0xff;
			uchar g = (curColor & 0xff00) >> 8;
			uchar r = (curColor & 0xff0000) >> 16;
			imRef(output, x, y).r = r;
			imRef(output, x, y).g = g;
			imRef(output, x, y).b = b;
		}
	}
	tempBMP->UnlockBits(&bitmapData);
	delete tempBMP;
	if (s_Token == 0UL) {
		Gdiplus::GdiplusShutdown(gdiplusToken);
	}
	return output;
}

image<rgb>* LoadImageGenericRGB(const char* filename,imgFormat* format)
{
	image<rgb>* output=0;
	if(format)*format=INVALID;
	output=loadPPMNoThrow(filename);
	if(output)
	{
		if(format)*format=PNM;
		return output;
	}
	image<uchar>*grey=loadPGMNoThrow(filename);
	if(grey)
	{
		output=new image<rgb>(grey->width(),grey->height(),false);
		for(int x=0;x<grey->width();x++)
		{
			for(int y=0;y<grey->height();y++)
			{
				imRef(output,x,y).r=imRef(grey,x,y);
				imRef(output,x,y).g=imRef(grey,x,y);
				imRef(output,x,y).b=imRef(grey,x,y);
			}
		}
		delete grey;
		if(format)*format=PNM;
		return output;
	}
	return LoadGDIPlusToRGB(filename, format);
}

image<double>* LoadImageGenericMono(const char* filename,imgFormat* format)
{
	image<double>* output=0;
	if(format)*format=INVALID;
	img::image<img::rgb>*colour=loadPPMNoThrow(filename);
	if(colour)
	{
		output=new image<double>(colour->width(),colour->height(),false);
		for(int x=0;x<colour->width();x++)
		{
			for(int y=0;y<colour->height();y++)
			{
				imRef(output,x,y)=YFromRGB(imRef(colour,x,y));
			}
		}
		delete colour;
		if(format)*format=PNM;
		return output;
	}
	image<uchar>*grey=loadPGMNoThrow(filename);
	if(grey)
	{
		output=new image<double>(grey->width(),grey->height(),false);
		for(int x=0;x<grey->width();x++)
		{
			for(int y=0;y<grey->height();y++)
			{
				imRef(output,x,y)=(double)imRef(grey,x,y);
			}
		}
		delete grey;
		if(format)*format=PNM;
		return output;
	}
	colour = LoadGDIPlusToRGB(filename, format);
	output = new image<double>(colour->width(), colour->height(), false);
	for (int x = 0; x < colour->width(); x++)
	{
		for (int y = 0; y < colour->height(); y++)
		{
			imRef(output, x, y) = YFromRGB(imRef(colour, x, y));
		}
	}
	delete colour;
	return output;
}

void SaveImageGeneric(const image<rgb>* img,const char* filename,imgFormat format, int jpegQuality)
{
	const WCHAR* mimetype=0;
	switch(format)
	{
		case PNM:
			{
				savePPM(img,filename);
			}
			return;
		case BMP:
			{
				mimetype=L"image/bmp";
			}
			break;
		case GIF:
			{
				mimetype=L"image/gif";
			}
			break;
		case JPEG:
			{
				mimetype=L"image/jpeg";
			}
			break;
		case PNG:
			{
				mimetype=L"image/png";
			}
			break;
		case TIFF:
			{
				mimetype=L"image/tiff";
			}
			break;
		case INVALID:
		default:
			return;
	}
	Gdiplus::GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR gdiplusToken;
	Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);
	size_t width=img->width();
	size_t height=img->height();
	Gdiplus::Bitmap*  bmp = new Gdiplus::Bitmap(static_cast<int>(width),static_cast<int>(height),PixelFormat24bppRGB);
	for(size_t x=0;x<width;x++)
	{
		for(size_t y=0;y<height;y++)
		{
			Gdiplus::Color color((BYTE)imRef(img,x,y).r,(BYTE)imRef(img,x,y).g,(BYTE)imRef(img,x,y).b);
			bmp->SetPixel(static_cast<int>(x),static_cast<int>(y),color);
		}
	}

	CLSID   encoderClsid;
	GetEncoderClsid(mimetype, &encoderClsid);

	int nLen = MultiByteToWideChar(CP_ACP,0,filename,-1,NULL,NULL);
	LPWSTR lpszW = new WCHAR[nLen];
	MultiByteToWideChar(CP_ACP,0,filename,-1,lpszW,nLen);
	if (format == imgFormat::JPEG) {
		Gdiplus::EncoderParameters params;
		params.Count = 1;
		params.Parameter[0].Guid = Gdiplus::EncoderQuality;
		params.Parameter[0].Type = Gdiplus::EncoderParameterValueTypeLong;
		params.Parameter[0].NumberOfValues = 1;
		long quality = static_cast<long>(jpegQuality);
		params.Parameter[0].Value = &quality;
		Gdiplus::Status stat = bmp->Save(lpszW, &encoderClsid, &params); // note default JPEG quality is 75
	} else {
		Gdiplus::Status stat = bmp->Save(lpszW, &encoderClsid, NULL);
	}
	delete[] lpszW;
	delete bmp;

	Gdiplus::GdiplusShutdown(gdiplusToken);
}

void SaveImageGeneric(const image<double>* img,const char* filename,imgFormat format)
{
	size_t width=img->width();
	size_t height=img->height();
	image<rgb>* tmpImg=new image<rgb>(width,height,false);
	double min,max;
	min=max=imRef(img,0,0);
	for(size_t x=0;x<width;x++)
	{
		for(size_t y=0;y<height;y++)
		{
			double value=imRef(img,x,y);
			if(value<min)min=value;
			if(value>max)max=value;
		}
	}
	for(size_t x=0;x<width;x++)
	{
		for(size_t y=0;y<height;y++)
		{
			uchar value=(uchar)std::clamp(255.0*((imRef(img,x,y)-min)/(max-min)),0.0,255.0);
			imRef(tmpImg,x,y).r=value;
			imRef(tmpImg,x,y).g=value;
			imRef(tmpImg,x,y).b=value;
		}
	}
	SaveImageGeneric(tmpImg,filename,format);
	delete tmpImg;
}

}
