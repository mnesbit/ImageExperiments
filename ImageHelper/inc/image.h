#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "framework.h"

#include <cstring>

namespace img{

typedef unsigned char uchar;

typedef struct {
	uchar r;
	uchar g;
	uchar b;
} rgb;

inline bool operator==(const rgb &a, const rgb &b) {
  return ((a.r == b.r) && (a.g == b.g) && (a.b == b.b));
}

inline rgb operator-(const rgb &a, const rgb &b) {
	rgb retval{};
	retval.r=(uchar)((int)a.r-(int)b.r+127);
	retval.g=(uchar)((int)a.g-(int)b.g+127);
	retval.b=(uchar)((int)a.b-(int)b.b+127);
	return retval;
}

typedef struct { double y, u, v; } yuv;
typedef struct { double X,Y,Z; } XYZ;

template <class T>
class image {
 public:
  /* create an image */
  image(const size_t width, const size_t height, const bool init = true);

  /* delete an image */
  ~image();

  /* init an image */
  void init(const T &val);

  /* copy an image */
  image<T> *copy() const;

  /* get the width of an image. */
  size_t width() const { return w; }

  /* get the height of an image. */
  size_t height() const { return h; }

  /* image data. */
  T *data;

  /* row pointers. */
  T **access;

 private:
  size_t w, h;
};

template <class T>
image<T>::image(const size_t width, const size_t height, const bool init) {
  w = width;
  h = height;
  data = new T[w * h];  // allocate space for image data
  access = new T*[h];   // allocate space for row pointers

  // initialize row pointers
  for (int i = 0; i < h; i++)
    access[i] = data + (i * w);

  if (init)
    memset(data, 0, w * h * sizeof(T));
}

template <class T>
image<T>::~image() {
  delete [] data;
  delete [] access;
}

template <class T>
image<T> *image<T>::copy() const {
  image<T> *im = new image<T>(w, h, false);
  memcpy(im->data, data, w * h * sizeof(T));
  return im;
}

/* use imRef to access image data. */
#define imRef(im, x, y) (im->access[y][x])

/* use imPtr to get pointer to image data. */
#define imPtr(im, x, y) &(im->access[y][x])

template <class T>
void image<T>::init(const T &val) {
  T *ptr = imPtr(this, 0, 0);
  T *end = imPtr(this, w-1, h-1);
  while (ptr <= end)
    *ptr++ = val;
}

}

#endif //_IMAGE_H_
