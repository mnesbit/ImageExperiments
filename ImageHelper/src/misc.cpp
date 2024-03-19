#include "../inc/misc.h"

namespace img
{

yuv YUVFromRGB(rgb in)
{
	yuv retval;
	double Y=(RED_WEIGHT*(double)in.r+GREEN_WEIGHT*(double)in.g+BLUE_WEIGHT*(double)in.b);
	retval.y=Y;
	retval.u=((U_WEIGHT/(1.0-BLUE_WEIGHT))* ((double)in.b-Y));
	retval.v=((V_WEIGHT/(1.0-RED_WEIGHT))*((double)in.r-Y));
	return retval;
}

rgb RGBFromYUV(yuv in)
{
	rgb retval;
	retval.r=(uchar)(in.y+1.13983*in.v);
	retval.g=(uchar)(in.y- 0.39466*in.u-0.58060*in.v);
	retval.b=(uchar)(in.y+ 2.03211*in.u);
	return retval;
}

double YFromRGB(rgb in)
{
	return (RED_WEIGHT*(double)in.r+GREEN_WEIGHT*(double)in.g+BLUE_WEIGHT*(double)in.b);
}

}