#include "render.h"
#include "image.h"
#include "parallel.h"
#include <vector>
#include <string>
#include <thread>

int main(int argc, char *argv[]) {
    std::vector<std::string> parameters;
    std::string outputname;
    int num_threads = std::thread::hardware_concurrency();
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-t") {
            num_threads = std::stoi(std::string(argv[++i]));
        } else if (std::string(argv[i]) == "-o") {
            outputname = std::string(argv[++i]);
        } else {
            parameters.push_back(std::string(argv[i]));
        }
    }

    parallel_init(num_threads);

    Image3 img = render_img(parameters);
    imwrite(outputname, img);

    parallel_cleanup();

    return 0;
}
