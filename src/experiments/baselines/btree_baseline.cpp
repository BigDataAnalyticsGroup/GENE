#include "util/Dataset.hpp"
#include "util/Workload.hpp"
#include "hp/Partition.hpp"
#include "hp/Index.hpp"

#include <iostream>
#include <fstream>
#include <chrono>
#include <ratio>
#include <cmath>
#include <numeric>

using key_type = uint64_t;
using mapped_type = uint64_t;

std::vector<long double> testIndex(HorizontalPartition<key_type, mapped_type>* hp, Dataset<key_type, mapped_type>* ds, Workload<key_type, mapped_type>* wl, unsigned int repetitions) {
    ds->sort();
    for (auto it=ds->begin(); it!=ds->end(); ++it) {
        try {
            hp->insert(it->first, it->second);
        } catch(CapacityException& e) {}
        catch(DuplicateKeyException& e) {}
    }

    std::vector<long double> durations;
    for (unsigned int r=0; r<repetitions; r++) {
        uint64_t result = 0;
        auto start = std::chrono::high_resolution_clock::now();
        for (auto it=wl->begin(); it!=wl->end(); ++it) {
            switch((*it)->workload_kind()) {
                case Workload_kind::Point:
                    {
                        auto p = static_cast<PointQuery<key_type, mapped_type>*>(*it);
                        result += *(hp->pointQ(p->key()));
                        break;
                    }
                case Workload_kind::Range:
                    {
                        auto r = static_cast<RangeQuery<key_type, mapped_type>*>(*it);
                        hp->rangeQ(result, r->lower(), r->upper());
                        break;
                    }
                case Workload_kind::Insert:
                    try {
                        auto i = static_cast<Insertion<key_type, mapped_type>*>(*it);
                        hp->insert(i->key(), i->value());
                    } catch(CapacityException& e) {}
                    break;
                default:
                    throw std::invalid_argument("Unknown Workload_type");
                }
        }
        auto end = std::chrono::high_resolution_clock::now();
        durations.push_back((std::chrono::duration<long double, std::milli>(end - start)).count());
    }
    return durations;
}

int main(int argc, char** argv) {
    if (argc < 5 || argc > 7)
        throw std::invalid_argument("Usage: btree_baseline <dataset> <workload> <fanout> <leafsize> <repetitions=50> <prefix=>");
    std::string filepath_ds = argv[1];
    std::string filepath_wl = argv[2];
    unsigned int fanout = std::stoi(argv[3]);
    unsigned int leafsize = std::stoi(argv[4]);

    auto ds = Dataset<key_type, mapped_type>::fromFile(filepath_ds);
    auto wl = Workload<key_type, mapped_type>::from_file(ds, filepath_wl);

    unsigned int repetitions = 50;
    if (argc == 6)
        repetitions = std::stoi(argv[5]);
    std::string prefix = "";
    if (argc == 7)
	    prefix = argv[6];

    std::cout << "Btree baseline" << std::endl;
    std::cout << "Dataset size: " << ds->size() << std::endl;
    std::cout << "Fanout: " << fanout << std::endl;
    std::cout << "Leafsize: " << leafsize << std::endl;
    std::cout << "Measurement repetitions: " << repetitions << std::endl;
    std::cout << "Prefix: " << prefix << std::endl;
    
    std::ofstream csv;
    if (prefix.empty())
    	csv.open("btree_baseline.csv");
    else
	    csv.open(prefix+"baseline.csv");
    csv << "Run" << "," << "Duration" << std::endl;
    
    auto btree_ptr = bulkload_btree<key_type, mapped_type>(ds->begin(), ds->end(), fanout, leafsize);
    HorizontalPartition<key_type, mapped_type>* btree = btree_ptr.release();
    std::ofstream f;
    if (prefix.empty())
    	f.open("btree_baseline.dot");
    else
	    f.open(prefix + "baseline.dot");
    btree->to_graphviz(f, false, false);
    f.close();

    auto durations = testIndex(btree, ds, wl, repetitions);
    for (unsigned int j=0; j<durations.size(); j++)
        csv << j << "," << durations[j] << std::endl;
    std::sort(durations.begin(), durations.end());
    long double medianDuration = repetitions % 2 == 0 ? (durations[repetitions / 2 - 1] + durations[repetitions / 2]) / 2 : durations[repetitions / 2];
    std::cout << "Median duration: " << medianDuration << std::endl;
    long double sum = std::accumulate(durations.begin(), durations.end(), 0.0);
    long double meanDuration = sum / durations.size();
    long double squared_sum = std::inner_product(durations.begin(), durations.end(), durations.begin(), 0.0);
    long double devDuration = std::sqrt(squared_sum / durations.size() - meanDuration * meanDuration);
    std::cout << "Mean duration: " << meanDuration << std::endl;
    std::cout << "Standard deviation: " << devDuration << std::endl;

    delete(btree);
    delete(ds);
    delete(wl);
    csv.close();
}
