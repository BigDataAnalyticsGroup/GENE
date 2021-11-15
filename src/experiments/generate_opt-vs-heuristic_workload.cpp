#include <cmath>
#include <err.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <util/Dataset.hpp>
#include <util/Workload.hpp>

using key_type = uint64_t;
using mapped_type = uint64_t;
using size_type = std::size_t;

int main(int argc, char* argv[]) {
    if (argc != 3 && argc != 4)
        err(EXIT_FAILURE, "Usage: '%s <path/to/dataset> <size> <path/to/workload/file>'", argv[0]);

    std::string path_to_dataset = argv[1];
    size_type num_queries = std::atoi(argv[2]);
    std::string path_to_wkl_file = argv[3];

    size_type num_point_p1 = 0.2 * num_queries;
    size_type num_point_p2 = 0.1 * num_queries;
    size_type num_range_p2 = 0.2 * num_queries;
    size_type num_point_p3 = 0.5 * num_queries;

    Dataset<key_type, mapped_type> *dataset = Dataset<key_type, mapped_type>::fromFile(path_to_dataset);
    Workload<key_type, mapped_type> workload(dataset);

    size_type n = dataset->size();
    size_type partition_first_border = 0.1 * n;
    size_type partition_second_border = 0.85 * n;

    /* First partition. */
    workload.add_point_uniform_index_based(num_point_p1, {0, partition_first_border - 1});
    /* Second partition. */
    workload.add_point_uniform_index_based(num_point_p2, {partition_first_border, partition_second_border - 1});
    workload.add_range_uniform_index_based(num_range_p2, 0.00001, {partition_first_border, partition_second_border -  1});
    /* Third partition. */
    workload.add_point_uniform_index_based(num_point_p3, {partition_second_border, n - 1});

    workload.to_file(path_to_wkl_file);

    delete(dataset);
    return 0;
}


