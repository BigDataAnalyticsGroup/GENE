#include "art/art_primary_lb.h"
#include "cxxopts.hpp"
#include "optimized/Index.hpp"
#include "pgm/pgm_index.hpp"
#include "util/Dataset.hpp"
#include "util/Workload.hpp"
#include <cfenv>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>
#include <tlx/container.hpp>

#define NO_SCAN

using namespace tlx;
using namespace std::chrono;

using size_type = std::size_t;
using key_type = uint64_t;
using mapped_type = uint64_t;
using value_type = std::pair<key_type, mapped_type>;
using hp_type = opt::HorizontalPartition<key_type, mapped_type>;

/* Global result variable to prevent optimizing away lookups. */
std::size_t result = 0;

/*======================================================================================================================
 * Configuration parameters
 *====================================================================================================================*/
const size_type NUM_RUNS = 5;
const size_type inner_slot_size = 16;
const size_type leaf_slot_size = 32;

/*======================================================================================================================
 * TLX B+Tree settings
 *====================================================================================================================*/
struct traits : tlx::btree_default_traits<key_type, mapped_type> {
    static const bool self_verify = false;
    static const bool debug = false;
    static const int inner_slots = inner_slot_size;
    static const int leaf_slots = leaf_slot_size;
    //static const size_t binsearch_threshold = 256;
};

/*======================================================================================================================
 * Main
 *====================================================================================================================*/
int main(int argc, char *argv[]) {
    fesetround(FE_DOWNWARD); // required by PGM to produce correct results


// TODO: exchange with cxxopts
    /* Argparser */
    cxxopts::Options options("optimized_vs_heuristic", "experiment comparing state-of-the-art indexes to GENE");

    /* Add program arguments. */
    options.add_options()
        ("dataset", "path to dataset", cxxopts::value<std::string>())
        ("workload", "path to workload", cxxopts::value<std::string>())
        ("output", "path to output", cxxopts::value<std::string>())
        ("h,help", "Usage")
    ;

    auto cl_parameter = options.parse(argc, argv);

    if (cl_parameter.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }

    /* Define the input values. */
    const std::string dataset_path = cl_parameter["dataset"].as<std::string>();
    const std::string workload_path = cl_parameter["workload"].as<std::string>();
    const std::string output_path = cl_parameter["output"].as<std::string>();

    /* Random Generator. */
    std::mt19937_64 g(42);

    /* Output file. */
    std::ofstream out(output_path, std::ios_base::app);
    if (not out) {
        std::cerr << "Failed to open the output file '" << output_path << "'" << std::endl;
        return (EXIT_FAILURE);
    }

    /* Data. */
    Dataset<key_type, mapped_type> *dataset = Dataset<key_type, mapped_type>::fromFile(dataset_path);

    /* Workload. */
    Workload<key_type, mapped_type> *workload = Workload<key_type, mapped_type>::from_file(dataset, workload_path);

    /* Define different variants of index structures.
     * 1) TLX B+Tree
     * 2) ARTful
     * 3) PGM
     * 4) GENE
     */
    /* Create 1) TLX B+Tree. */
    btree_map<key_type, mapped_type, std::less<key_type>, traits> tlx_btree;
    /* Note: When bulkloading a `tlx::btree_map` the data needs to be *duplicate free*. */
    tlx_btree.bulk_load(dataset->begin(), dataset->end());

    /* Create 2) ARTful. */
    ARTPrimaryLB<1> art;
    art.Build(dataset->entries());

    /* Create 3) PGM. */
    std::vector<key_type> data_pgm;
    for (auto d : *dataset) { data_pgm.push_back(d.first); }
    const int epsilon = 128; // space-time trade-off parameter
    pgm::PGMIndex<key_type, epsilon> pgm_index(data_pgm);

    /* Create 4) GENE. */
    auto gene = opt::Index<key_type, mapped_type>::Bulkload(dataset->begin(), dataset->end(),
                                                            inner_slot_size + 1, leaf_slot_size);

#define OUTPUT(N, DIST, WORKLOAD_TYPE, WORKLOAD_SIZE, INDEX_TYPE, INNER_SIZE, LEAF_SIZE, NUM_RUN, TIME) { \
    out << N << "," \
        << DIST << "," \
        << WORKLOAD_TYPE << "," \
        << WORKLOAD_SIZE << "," \
        << INDEX_TYPE << "," \
        << INNER_SIZE << "," \
        << LEAF_SIZE << "," \
        << NUM_RUN << "," \
        << TIME << std::endl; \
}

    /* Compute checksum for given workload. */
    uint64_t checksum = 0;
    for (auto it = workload->begin(), end = workload->end(); it != end; ++it) {
        switch((*it)->workload_kind()) {
            case Workload_kind::Point:
                {
                    auto p = static_cast<PointQuery<key_type, mapped_type>*>(*it);
                    auto key = p->key();
                    auto it = std::lower_bound(dataset->begin(), dataset->end(), key, [](auto const &lhs, auto const &k) { return lhs.first < k; });
                    if (it != dataset->end() and it->first == key) checksum += it->second;
                    break;
                }
            case Workload_kind::Range:
                {
                    auto r = static_cast<RangeQuery<key_type, mapped_type>*>(*it);
                    auto lower = r->lower();
                    auto it = std::lower_bound(dataset->begin(), dataset->end(), lower, [](auto const &lhs, auto const &k) { return lhs.first < k; });
                    if (it != dataset->end()) checksum += it->second;
                    break;
                }
            default:
                    throw std::invalid_argument("Unknown Workload_type");
        }
    }
    std::cout << "Checksum: " << checksum << std::endl;

    /* 1) TLX B+Tree. */
    for (auto i = 0; i != NUM_RUNS; ++i) {
        workload->shuffle();
        result = 0;

        auto start = steady_clock::now();
        for (auto it = workload->begin(), end = workload->end(); it != end; ++it) {
            switch((*it)->workload_kind()) {
                case Workload_kind::Point:
                    {
                        auto p = static_cast<PointQuery<key_type, mapped_type>*>(*it);
                        auto it = tlx_btree.find(p->key());
                        if (it != tlx_btree.end()) result += it->second;
                        break;
                    }
                case Workload_kind::Range:
                    {
                        auto r = static_cast<RangeQuery<key_type, mapped_type>*>(*it);
                        auto upper = r->upper();
                        auto data_size = dataset->size();
                        auto lb_it = tlx_btree.lower_bound(r->lower());
                        if (lb_it == tlx_btree.end()) continue;
#ifndef NO_SCAN
                        auto idx = lb_it->second; // offset into data array
                        while (idx != data_size && (*dataset)[idx].first <= upper) {
                            ++result;
                            ++idx;
                        }
#else
                        result += lb_it->second;
#endif
                        break;
                    }
                case Workload_kind::Insert:
                    std::cerr << "insert not supposed to be in this workload" << std::endl;
                    break;
                default:
                    throw std::invalid_argument("Unknown Workload_type");
            }
        }
        auto stop = steady_clock::now();
        auto time = duration_cast<nanoseconds>(stop - start).count();
        /* Validate checksum. */
        std::cout << "TLX Checksum: " << result << std::endl;
        assert(result == checksum);
        OUTPUT(dataset->size(), toString(dataset->distribution()), "poc", workload->size(), "TLX B-tree", traits::inner_slots, traits::leaf_slots, i, time);
    }

    /* 2) ARTful. */
    for (auto i = 0; i != NUM_RUNS; ++i) {
        workload->shuffle();
        result = 0;

        auto start = steady_clock::now();
        for (auto it = workload->begin(), end = workload->end(); it != end; ++it) {
            switch((*it)->workload_kind()) {
                case Workload_kind::Point:
                    {
                        auto p = static_cast<PointQuery<key_type, mapped_type>*>(*it);
                        auto v = art.find(p->key());
                        if (v) result += *v;
                        break;
                    }
                case Workload_kind::Range:
                    {
                        auto r = static_cast<RangeQuery<key_type, mapped_type>*>(*it);
                        auto upper = r->upper();
                        auto data_size = dataset->size();
                        auto lb_opt = art.lower_bound(r->lower());
                        if (not lb_opt) continue;
#ifndef NO_SCAN
                        auto idx = *lb_opt;
                        while (idx != data_size && (*dataset)[idx].first <= upper) {
                            ++result;
                            ++idx;
                        }
#else
                        result += *lb_opt;
#endif
                        break;
                    }
                case Workload_kind::Insert:
                    std::cerr << "insert not supposed to be in this workload" << std::endl;
                    break;
                default:
                    throw std::invalid_argument("Unknown Workload_type");
                }
        }
        auto stop = steady_clock::now();
        auto time = duration_cast<nanoseconds>(stop - start).count();
        /* Validate checksum. */
        std::cout << "ART Checksum: " << result << std::endl;
        assert(result == checksum);
        OUTPUT(dataset->size(), toString(dataset->distribution()), "poc", workload->size(), "ART", 0, 0, i, time);
    }

    /* 3) PGM. */
    for (auto i = 0; i != NUM_RUNS; ++i) {
        workload->shuffle();
        result = 0;

        auto start = steady_clock::now();
        for (auto it = workload->begin(), end = workload->end(); it != end; ++it) {
            switch((*it)->workload_kind()) {
                case Workload_kind::Point:
                    {
                        auto p = static_cast<PointQuery<key_type, mapped_type>*>(*it);
                        auto key = p->key();
                        auto range = pgm_index.search(key);
                        auto lo = data_pgm.begin() + range.lo;
                        auto hi = data_pgm.begin() + range.hi;
                        auto pos = std::lower_bound(lo, hi, key);
                        if (*pos == key) result += std::distance(data_pgm.begin(), pos);
                        break;
                    }
                case Workload_kind::Range:
                    {
                        auto r = static_cast<RangeQuery<key_type, mapped_type>*>(*it);
                        auto lower = r->lower();
                        auto upper = r->upper();
                        auto range = pgm_index.search(lower);
                        auto lo = data_pgm.begin() + range.lo;
                        auto hi = data_pgm.begin() + range.hi;
                        auto pos = std::lower_bound(lo, hi, lower);
#ifndef NO_SCAN
                        auto idx = std::distance(data_pgm.begin(), lb_it);
                        while (idx != data_size && data_pgm[idx] <= upper) {
                            ++result;
                            ++idx;
                        }
#else
                        result += std::distance(data_pgm.begin(), pos);;
                        break;
#endif
                    }
                case Workload_kind::Insert:
                    std::cerr << "insert not supposed to be in this workload" << std::endl;
                    break;
                default:
                    throw std::invalid_argument("Unknown Workload_type");
                }
        }
        auto stop = steady_clock::now();
        auto time = duration_cast<nanoseconds>(stop - start).count();
        /* Validate checksum. */
        std::cout << "PGM Checksum: " << result << std::endl;
        assert(result == checksum);
        OUTPUT(dataset->size(), toString(dataset->distribution()), "poc", workload->size(), "PGM", 0, 0, i, time);
    }

    /* 4) GENE. */
    for (auto i = 0; i != NUM_RUNS; ++i) {
        workload->shuffle();
        result = 0;

        auto start = steady_clock::now();
        for (auto it = workload->begin(), end = workload->end(); it != end; ++it) {
            switch((*it)->workload_kind()) {
                case Workload_kind::Point:
                    {
                        auto p = static_cast<PointQuery<key_type, mapped_type>*>(*it);
                        auto v = gene.find(p->key());
                        if (v) result += *v;
                        break;
                    }
                case Workload_kind::Range:
                    {
                        auto r = static_cast<RangeQuery<key_type, mapped_type>*>(*it);
                        auto upper = r->upper();
                        auto data_size = dataset->size();
                        auto lb_opt = gene.lower_bound(r->lower());
                        if (not lb_opt) continue;
#ifndef NO_SCAN
                        auto idx = *lb_it;
                        while (idx != data_size && (*dataset)[idx].first <= upper) {
                            ++result;
                            ++idx;
                        }
#else
                        result += *lb_opt;
#endif
                        break;
                    }
                case Workload_kind::Insert:
                    std::cerr << "insert not supposed to be in this workload" << std::endl;
                    break;
                default:
                    throw std::invalid_argument("Unknown Workload_type");
                }
        }
        auto stop = steady_clock::now();
        auto time = duration_cast<nanoseconds>(stop - start).count();
        /* Validate checksum. */
        std::cout << "GENE Checksum: " << result << std::endl;
        assert(result == checksum);
        OUTPUT(dataset->size(), toString(dataset->distribution()), "poc", workload->size(), "GENE", 3, 0, i, time);
    }

    /* Delete data & workload. */
    delete(dataset);
    delete(workload);

    /* Close output files. */
    out.close();

    return 0;
}
#undef OUTPUT
