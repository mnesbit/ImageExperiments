#include "../inc/misc.h"

namespace img
{

	yuv YUVFromRGB(rgb in)
	{
		return YUVFromRGB((double)in.r, (double)in.g, (double)in.b);
	}

	yuv YUVFromRGB(double red, double green, double blue)
	{
		double Y = (RED_WEIGHT * red + GREEN_WEIGHT * green + BLUE_WEIGHT * blue);
		yuv retval{
			.y = Y,
			.u = ((U_WEIGHT / (1.0 - BLUE_WEIGHT)) * (blue - Y)),
			.v = ((V_WEIGHT / (1.0 - RED_WEIGHT)) * (red - Y))
		};
		return retval;
	}

	rgb RGBFromYUV(yuv in)
	{
		rgb retval{
			.r = (uchar)(in.y + 1.13983 * in.v),
			.g = (uchar)(in.y - 0.39466 * in.u - 0.58060 * in.v),
			.b = (uchar)(in.y + 2.03211 * in.u)
		};
		return retval;
	}

	double YFromRGB(rgb in)
	{
		return (RED_WEIGHT * (double)in.r + GREEN_WEIGHT * (double)in.g + BLUE_WEIGHT * (double)in.b);
	}

}