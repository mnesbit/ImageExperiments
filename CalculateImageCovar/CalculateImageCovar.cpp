#include <iostream>
#include <filesystem>
#include <random>
#include "../ImageHelper/inc/imgloader.h"
#include "../SimpleMatrix/inc/covariance.h"
#include "../ImageHelper/inc/misc.h"

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
            }
            else if ((sb.st_mode & S_IFDIR) == S_IFDIR) {
                const auto subdirectory = readFiles(path);
                results.insert(results.end(), subdirectory.cbegin(), subdirectory.cend());
            }
        }
    }
    return results;
}

int main(int argc, char* argv[]) {
    if (argc != 5 && argc != 6) {
        std::cout << "usage: CalculateImageCovar.exe <patch size in pixels> <random patches per image> <image source folder> <output matrix file name> [<output PNG file name for graphical version>]" << std::endl;
        return -1;
    }
    int N = atoi(argv[1]);
    int patchesPerImage = atoi(argv[2]);
    std::string sourceFolder(argv[3]);
    std::string outputFile(argv[4]);
    std::cout << std::format("running with patch size={}x{} sampled {} times from each image in '{}' and output to '{}'", N, N, patchesPerImage, sourceFolder, outputFile) << std::endl;
    if (argc == 6) {
        std::cout << std::format("also writing covariance as greyscale image to {}", argv[5]) << std::endl;
    }
    Startup();
    const auto files = readFiles(sourceFolder);
    math::CovarianceCalculator covar;
    math::Vector patch(N * N);
    std::mt19937 rand;
    rand.seed(123456);
    for (size_t count = 0; count < files.size(); ++count) {
        std::string file = files[count];
        std::string lowerCaseName(file.size(), 0);
        transform(file.cbegin(), file.cend(), lowerCaseName.begin(), ::tolower);
        if (lowerCaseName.ends_with(".jpg")) {
            std::cout << std::format("processing: {} {}%", file, (100L * count)/files.size()) << std::endl;
            imgFormat format;
            image<rgb>* image = LoadImageGenericRGB(file.c_str(), &format);
            if (image->width() < N || image->height() < N) {
                delete image;
                continue;
            }
            for (int reps = 0; reps < patchesPerImage; ++reps) {
                int x = rand() % (image->width() - N);
                int y = rand() % (image->height() - N);
                for (int offx = 0; offx < N; ++offx) {
                    for (int offy = 0; offy < N; ++offy) {
                        patch[offx + N * offy] = YFromRGB(imRef(image, x + offx, y + offy));
                    }
                }
                covar.Update(patch);
            }
            delete image;
            covar.Covariance().Save(outputFile.c_str());
        }
    }
    if (argc == 6) {
        math::Matrix outputCovar = covar.Covariance();
        image<double>* img = new image<double>(N * N, N * N, false);
        for (int x = 0; x < N * N; ++x) {
            int x1 = x % N;
            int y1 = x / N;
            for (int y = 0; y < N * N; ++y) {
                int x2 = y % N;
                int y2 = y / N;
                imRef(img, x1 * N + x2, y1 * N + y2) = outputCovar[x][y];
            }
        }
        std::cout << "writing covariance image file to: " << argv[5] << std::endl;
        SaveImageGeneric(img, argv[5], imgFormat::PNG);
        delete img;
    }
    Shutdown();
    return 0;
}