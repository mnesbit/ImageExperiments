#include "../inc/misc.h"
#include <algorithm>

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
		return RGBFromYUV(in.y, in.u, in.v);
	}

	rgb RGBFromYUV(double y, double u, double v) {
		rgb retval{
			.r = static_cast<uchar>(std::clamp(std::round(y + 1.13983 * v), 0.0, 255.0)),
			.g = static_cast<uchar>(std::clamp(std::round(y - 0.39466 * u - 0.58060 * v), 0.0, 255.0)),
			.b = static_cast<uchar>(std::clamp(std::round(y + 2.03211 * u), 0.0, 255.0))
		};
		return retval;

	}

	double YFromRGB(rgb in)
	{
		return (RED_WEIGHT * (double)in.r + GREEN_WEIGHT * (double)in.g + BLUE_WEIGHT * (double)in.b);
	}

}