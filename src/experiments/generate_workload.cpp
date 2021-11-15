#include <cstdlib>
#include <err.h>
#include <stdlib.h>
#include <set>

#include "util/Dataset.hpp"
#include "util/Workload.hpp"


using key_type = uint64_t;
using size_type = std::size_t;


/* https://stackoverflow.com/questions/8520560/get-a-file-name-from-a-path */
std::vector<std::string> splitpath( const std::string& str , const std::set<char> delimiters)
{
    std::vector<std::string> result;

    char const* pch = str.c_str();
    char const* start = pch;

    for(; *pch; ++pch) {
        if (delimiters.find(*pch) != delimiters.end()) {
            if (start != pch) {
                std::string str(start, pch);
                result.push_back(str);
            }
            else {
                result.push_back("");
            }
            start = pch + 1;
        }
    }

    result.push_back(start);
    return result;
}

/**
 * Generate mixed workload for the given dataset.
 */
int main(int argc, char* argv[]) {
    if (argc != 5 && argc != 6)
        err(EXIT_FAILURE, "Usage: '%s <path/to/dataset> <size> <ratio> <sel> <path/to/result/file : optional>'\n\t- sel = [0,1]\n\t- ratio (percentage point queries) = [0,1]", argv[0]);
    std::string file = std::string(argv[1]);
    std::vector<std::string> path = splitpath(file, {'/'});
    std::string s = path.back();
    std::string dataset_name = s.substr(0, s.find(std::string(".")));
    size_type size = std::atoi(argv[2]);
    double ratio = std::atof(argv[3]);
    double sel = std::atof(argv[4]);
    uint64_t num_point = ratio * size;
    uint64_t num_range = size - num_point;

    Dataset<key_type, key_type> *data = Dataset<key_type, key_type>::fromFile(file);

    key_type dom_min = data->domain().min;
    key_type dom_max = data->domain().max;

    Workload<key_type, key_type> workload(data);
    workload.add_point_uniform(num_point, {dom_min, dom_max});
    workload.add_range_uniform_index_based(num_range, sel, {dom_min, dom_max});
    std::string out = "workloads/mix_" + dataset_name + "_" + std::to_string(size) + "_" + std::to_string(ratio) + "_" + std::to_string(sel) + ".wkl";
    if (argc == 6)
        out = std::string(argv[5]);
    workload.to_file(out);

    delete data;
}

