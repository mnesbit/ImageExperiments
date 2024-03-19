// ImageLoader.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <filesystem>
#include "../ImageHelper/inc/imgloader.h"
using namespace img;

std::vector<std::string> readFiles(std::string path) {
    struct stat sb;
    std::vector<std::string> results;

    for (const auto& entry : std::filesystem::directory_iterator(path)) {

        std::filesystem::path outfilename = entry.path();
        std::string outfilename_str = outfilename.string();
        const char* path = outfilename_str.c_str();

        if (stat(path, &sb) == 0) {
            if (!(sb.st_mode & S_IFDIR)) {
                results.push_back(outfilename_str);
            } else if ((sb.st_mode & S_IFDIR) == S_IFDIR) {
                const auto subdirectory = readFiles(path);
                results.insert(results.end(),subdirectory.cbegin(), subdirectory.cend());
            }
        }
    }
    return results;
}

int main()
{
    std::cout << "Hello World!\n";
	//imgFormat format;
	//image<rgb>* image = LoadImageGenericRGB("D:\\Test\\pool\\2007_001073.jpg", &format);
	////image<rgb>* image = LoadImageGenericRGB("C:\\Temp\\Barbara.pnm", &format);
	//SaveImageGeneric(image, "c:\\temp\\out.jpg", JPEG);
	//delete image;
    Startup();
    std::string path = "D:\\Treasure\\docs\\My Pictures\\";
    const auto files = readFiles(path);
    for (const auto& file : files) {
        std::string lowerCaseName(file.size(), 0);
        transform(file.cbegin(), file.cend(), lowerCaseName.begin(), ::tolower);
        if (lowerCaseName.ends_with(".jpg")) {
            std::cout << file << std::endl;
            imgFormat format;
            image<rgb>* image = LoadImageGenericRGB(file.c_str(), &format);
            delete image;
        }
    }
    Shutdown();
	return 0;

}
