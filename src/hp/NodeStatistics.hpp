#pragma once

#include "util/enum.hpp"

template<typename Key, typename Value>
struct NodeStatistics
{
    using size_type = std::size_t;
    using key_type = Key;
    using mapped_type = Value;

    public:

    /* ==========================================================================================
     * Constructor & Destructor
     * ========================================================================================*/

    NodeStatistics(PartitionType type, key_type key=key_type(), long double total_time=0, unsigned int num_point=0, 
            unsigned int num_range=0, unsigned int num_insert=0, unsigned int num_splits=0)
        : key(key),
          type(type),
          total_time(total_time),
          num_point(num_point),
          num_range(num_range),
          num_insert(num_insert),
          num_splits(num_splits)
    {}

    NodeStatistics(NodeStatistics* ns)
        : key(ns->key),
          type(ns->type),
          total_time(ns->total_time),
          num_point(ns->num_point),
          num_range(ns->num_range),
          num_insert(ns->num_insert),
          num_splits(ns->num_splits)
    {}

    /* ==========================================================================================
     * Manage executions of a node
     * ========================================================================================*/
    
    void add_execution(long double time, Workload_kind type) {
        total_time += time;
        switch(type) {
            case Workload_kind::Point:
                num_point++;
                break;
            case Workload_kind::Range:
                num_range++;
                break;
            case Workload_kind::Insert:
                num_insert++;
                break;
            default:
                std::cerr << "Unexpected workload type in NodeStatistics::add_execution";
                throw 666;
        }
    }

    void add_split() {
        num_splits++;
    }

    void add_penalty(long double time) {
        total_time += time;
    }

    void remove_execution(long double time, Workload_kind type) {
        total_time -= time;
        switch(type) {
            case Workload_kind::Point:
                num_point--;
                break;
            case Workload_kind::Range:
                num_range--;
                break;
            case Workload_kind::Insert:
                num_insert--;
                break;
            default:
                std::cerr << "Unexpected workload type in NodeStatistics::add_execution";
                throw 666;
        }
    }

    void reset_statistics() {
        total_time = 0;
        num_point = 0;
        num_range = 0;
        num_insert = 0;
    }

    void normalize_repetitions(uint64_t num_rep) {
        total_time /= num_rep;
        num_point  /= num_rep;
        num_range  /= num_rep;
        num_insert /= num_rep;
    }

    void add_results(const NodeStatistics<key_type, mapped_type>* other) {
        total_time += other->total_time;
        num_point += other->num_point;
        num_range += other->num_range;
        num_insert += other->num_insert;
        num_splits += other->num_splits;
    }

    long double approximate_costs() const {
        long double costs_point = 1;
        long double costs_range = 1;
        long double costs_insert = 1;
        long double costs_split = 1;

        switch (type) {
            case PartitionType::SoAPartition:
                costs_insert = 1;
                costs_split = 20 * costs_insert;
                break;
            case PartitionType::HashPartition:
                costs_insert = 3;
                costs_split = 15 * costs_insert;
                break;
            case PartitionType::TreePartition:
                costs_insert = 3;
                costs_split = 15 * costs_insert;
                break;
            default:
                break;
        }
        return num_point * costs_point + num_range * costs_range + num_insert * costs_insert + num_splits * costs_split;
    }


    /* ==========================================================================================
     * Getter & setter
     * ========================================================================================*/

    void set_key(key_type k) { key = k; }
    key_type get_key() const { return key; }
    void set_partition_type(PartitionType t) { type = t; }
    PartitionType get_partition_type() { return type; }
    long double get_total_time() const { 
        return approximate_costs();
    }
    unsigned int get_num_calls() const { return num_point + num_range + num_insert; }
    long double get_average_time() const { return total_time / get_num_calls(); }
    unsigned int get_num_point() const { return num_point; }
    unsigned int get_num_range() const { return num_range; }
    unsigned int get_num_insert() const { return num_insert; }
    unsigned int get_num_splits() const { return num_splits; }

    /* ==========================================================================================
     * Class methods
     * ========================================================================================*/
    void print(std::ostream& os=std::cout) const {
        os << "NodeStatistics(Key: " << get_key() << ", Total Time: " << get_total_time() << ", Total Calls: " << get_num_calls() <<  ", Point Queries: " << num_point << ", Range Queries: " << num_range << ", Inserts: " << num_insert << ", Splits: " << num_splits << ")";
    }

    friend std::ostream& operator<<(std::ostream& out, const NodeStatistics& ns) {
        ns.print(out);
        return out;
    }

    /* ==========================================================================================
     * Class members
     * ========================================================================================*/
    
    private:
        key_type key;
        PartitionType type;
        long double total_time;
        unsigned int num_point;
        unsigned int num_range;
        unsigned int num_insert;
        unsigned int num_splits;
};
