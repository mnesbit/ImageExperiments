#ifndef _PNM_FILE_H_
#define _PNM_FILE_H_

#include "framework.h"

#include <exception>
#include <memory>
#include "image.h"

namespace pnm {

class pnm_error:public std::exception { };

std::unique_ptr<img::image<img::uchar> > loadPBM(const char *name);
std::unique_ptr<img::image<img::uchar> > loadPBMNoThrow(const char *name);
void savePBM(const img::image<img::uchar> *im, const char *name);
std::unique_ptr<img::image<img::uchar> > loadPGM(const char *name) ;
std::unique_ptr<img::image<img::uchar> > loadPGMNoThrow(const char *name) ;
void savePGM(const img::image<img::uchar> *im, const char *name);
std::unique_ptr<img::image<img::rgb> > loadPPM(const char *name);
std::unique_ptr<img::image<img::rgb> > loadPPMNoThrow(const char *name);
void savePPM(const img::image<img::rgb> *im, const char *name);

template <class T>
void load_image(img::image<T> **im, const char *name);

template <class T>
void save_image(const img::image<T> *im, const char *name);

}

#endif //_PNM_FILE_H_
