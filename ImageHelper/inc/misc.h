#ifndef _MISC_H_
#define _MISC_H_

#include <cmath>
#include "image.h"

#define	RED_WEIGHT		0.299
#define GREEN_WEIGHT	0.587
#define BLUE_WEIGHT		0.114
#define U_WEIGHT		0.436
#define V_WEIGHT		0.615

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

template<class T>
class ArrayDeleter
{
public:
    void operator () (T* d) const
    {
        delete [] d;
    }
};

template <class T>
inline T abs(const T &x) { return (x > 0 ? x : -x); };

template <class T>
inline int sign(const T &x) { return (x >= 0 ? 1 : -1); };

template <class T>
inline T square(const T &x) { return x*x; };

template <class T>
inline bool check_bound(const T &x, const T&min, const T &max) {
  return ((x < min) || (x > max));
}

namespace img
{

typedef struct Point_Tag
{
	int x;
	int y;
} Point;

yuv YUVFromRGB(double red, double green, double blue);
yuv YUVFromRGB(rgb in);

rgb RGBFromYUV(yuv in);

double YFromRGB(rgb in);

}

#endif //_MISC_H_
