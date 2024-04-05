#include <iostream>

#include "../ImageHelper/inc/image.h"
#include "../ImageHelper/inc/imgloader.h"
#include "cluster.h"

int main(int argc, char* argv[])
{
    int rounds = atoi(argv[1]);
    img::image<img::rgb>* imgIn = img::LoadImageGenericRGB(argv[2]);
    cluster::Clusters cluster(imgIn);
    cluster.Process(abs(rounds));
    delete imgIn;
    img::image<img::rgb>* imgOut;
    if (rounds < 0) {
        imgOut  = cluster.Borders();
    } else {
        imgOut = cluster.Averages();
    }
    img::SaveImageGeneric(imgOut, argv[3], img::PNG);
    delete imgOut;
}
