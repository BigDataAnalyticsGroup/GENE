#include <cstring>
#include <err.h>
#include <stdlib.h>

#include "util/Dataset.hpp"


using key_type = uint64_t;
using size_type = std::size_t;

int main(int argc, char* argv[]) {
    if (argc != 3 && argc != 4)
        err(EXIT_FAILURE, "Usage: '%s <dataset> <size> <path/to/result/file : optional>'\n\t- dataset = [uni-dense, uni-sparse, books]", argv[0]);

    /* Output file. */
    std::string out_file = "data/" + std::string(argv[1]) + '_' + std::string(argv[2]) + ".data";
    if (argc == 4)
        out_file = std::string(argv[3]);

    size_type size = atoi(argv[2]);

    if (strcmp(argv[1], "uni-dense") == 0) {
        Dataset<key_type, key_type> * data = Dataset<key_type, key_type>::Create_Uniform_Dense(size, 0);
        data->toFile(out_file);
        delete(data);
    } else if (strcmp(argv[1], "uni-sparse") == 0) {
        Dataset<key_type, key_type> * data = Dataset<key_type, key_type>::Create_Uniform_Sparse(size, {0, std::numeric_limits<key_type>::max()});
        data->toFile(out_file);
        delete(data);
    } else if (strcmp(argv[1], "books") == 0) {
        Dataset<key_type, key_type> * data = Dataset<key_type, key_type>::Create_Books(size);
        data->toFile(out_file);
        delete(data);
    } else if (strcmp(argv[1], "osm") == 0) {
        Dataset<key_type, key_type> * data = Dataset<key_type, key_type>::Create_Osm(size);
        data->toFile(out_file);
        delete(data);
    } else {
        err(EXIT_FAILURE, "Wrong dataset provided!");
    }

    return 0;
}
