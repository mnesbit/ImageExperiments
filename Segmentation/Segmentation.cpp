#include <iostream>

#include "../ImageHelper/inc/image.h"
#include "../ImageHelper/inc/imgloader.h"
#include "cluster.h"

int main(int argc, char* argv[])
{
    if (argc != 4) {
        std::cout << "Segmentation.exe <cluster merge rounds> <input file> <output file>" << std::endl;
        return -1;
    }
    int rounds = atoi(argv[1]);
    std::unique_ptr<img::image<img::rgb> > imgIn = img::LoadImageGenericRGB(argv[2]);
    cluster::Clusters cluster(imgIn.get());
    cluster.Process(abs(rounds));
    imgIn.reset();
    std::unique_ptr<img::image<img::rgb> > imgOut;
    if (rounds < 0) {
        imgOut  = cluster.Borders();
    } else {
        imgOut = cluster.Averages();
    }
    img::SaveImageGeneric(imgOut.get(), argv[3], img::PNG);
    return 0;
}
