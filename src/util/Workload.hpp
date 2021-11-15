#pragma once

#include <algorithm>
#include <random>
#include <string>
#include <type_traits>
#include <cassert>
#include <fstream>
#include <sstream>

#include "util/Dataset.hpp"
#include "util/util.hpp"


#define IMPORT(WHAT) using WHAT = typename Base::WHAT

template<typename Key, typename T>
struct Workload_type
{
    using key_type = Key;
    using mapped_type = T;
    using size_type = std::size_t;

    protected:
    Workload_kind kind_;

    public:
    virtual ~Workload_type() {}

    Workload_kind workload_kind() const { return kind_; }
    virtual std::string to_string() const = 0;
};

template<typename Key, typename T>
struct PointQuery final : public Workload_type<Key, T>
{
    using Base = Workload_type<Key, T>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(size_type);

    private:
    const key_type key_;

    public:
    PointQuery(key_type key)
        : key_(key)
    {
        Base::kind_ = Workload_kind::Point;
    }

    const key_type key() const { return key_; }

    std::string to_string() const {
        std::stringstream ss;
        ss << "point:" << key();
        return ss.str();
    }
};

template<typename Key, typename T>
struct RangeQuery final : public Workload_type<Key, T>
{
    using Base = Workload_type<Key, T>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(size_type);

    private:
    const key_type lower_;
    const key_type upper_;

    public:
    RangeQuery(key_type lower, key_type upper)
        : lower_(lower)
        , upper_(upper)
    {
        Base::kind_ = Workload_kind::Range;
    }

    const key_type lower() const { return lower_; }
    const key_type upper() const { return upper_; }

    std::string to_string() const {
        std::stringstream ss;
        ss << "range:" << lower() << ":" << upper();
        return ss.str();
    }
};

template<typename Key, typename T>
struct Insertion final : public Workload_type<Key, T>
{
    using Base = Workload_type<Key, T>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(size_type);

    private:
    const key_type key_;
    const mapped_type value_;

    public:
    Insertion(key_type key, mapped_type value)
        : key_(key)
        , value_(value)
    {
        Base::kind_ = Workload_kind::Insert;
    }

    const key_type key() const { return key_; }
    const mapped_type value() const { return value_; }

    std::string to_string() const {
        std::stringstream ss;
        ss << "insert:" << key() << ":" << value();
        return ss.str();
    }
};

#undef IMPORT

template<typename Key, typename T>
struct Workload
{
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<const Key, T>;
    using size_type = std::size_t;
    using dataset_type = Dataset<key_type, mapped_type>;
    using wl_type = Workload_type<key_type, mapped_type>;
    using domain_range = typename dataset_type::domain_range;
    using index_range = typename util::Range<size_type>;

    private:
    std::vector<wl_type*> entries_;  //< workload entries
    dataset_type *dataset_;         //< underlying dataset
    uint64_t seed_;                 //< random generator seed
    std::mt19937_64 g_;             //< random generator
    bool has_inserts_ = false;       //< flag if workload contains inserts

    public:
    Workload(dataset_type *dataset, uint64_t seed = 21)
        : dataset_(dataset)
        , seed_(seed)
        , g_(std::mt19937_64(seed))
    {
        static_assert(std::is_integral<key_type>::value, "Integral type for Key required.");
        static_assert(std::is_integral<mapped_type>::value, "Integral type for T required.");
    }

    ~Workload() { for (auto e : entries_) delete e; }

    void add(wl_type *entry) {
        entries_.push_back(entry);
        has_inserts_ = has_inserts_ or entry->workload_kind() == Workload_kind::Insert;
    }

/* ============================================== POINT ==============================================================*/
    /** Create point queries for the underlying dataset and add them to the workload.
     *  The `n` elements are individually taken from `dataset_` starting at the very first index position. */
    void add_point_uniform_dense(size_type n) {
        assert(n <= dataset_->size());
        entries_.reserve(n);
        for (size_type i = 0; i < n; ++i)
            entries_.push_back(new PointQuery<key_type, mapped_type>((*dataset_)[i].first));
    }

    /** Create point queries for the underlying dataset and add them to the workload.
     *  The `n` elements are drawn uniformly at random in the key domain `dom` from `dataset_`.
     *  Note: the domain `dom` does refer to the key domain of the dataset. */
    void add_point_uniform(size_type n, domain_range dom) {
        entries_.reserve(n);
        std::uniform_int_distribution<key_type> distr(dom.min, dom.max);
        for (size_type i = 0; i < n; ++i)
            entries_.push_back(new PointQuery<key_type, mapped_type>(distr(g_)));
    }

    /** Create point queries for the underlying dataset which are definitely conainted, i.e. each point query hits
     *  exactly one element in the index, and add them to the workload.
     *  The `n` elements are drawn uniformly at random in the index domain `dom` from `dataset_`.
     *  Note: the index `range` does *not* refer to the key domain of the dataset but instead to the numbered elements in
     *  `dataset_`. For example, the dataset contains 10 elements, indexed from 0 to 9. The domain [1,3] restricts the
     *  point queries to be generated to the elements at position 1, 2, and 3 in the dataset. Therefore, the order of
     *  elements in `dataset_` matters! */
    void add_point_uniform_index_based(size_type n, index_range range) {
        entries_.reserve(n);
        assert(range.max < dataset_->size());
        std::uniform_int_distribution<size_type> distr(range.min, range.max);
        for (size_type i = 0; i < n; ++i)
            entries_.push_back(new PointQuery<key_type, mapped_type>((*dataset_)[distr(g_)].first));
    }

    /** Create point queries for the underlying dataset and add them to the workload.
     *  The `n` elements are drawn uniformly at random in the key domain from `dataset_`.
     *  Note: the domain `dom` does refer to the key domain of the dataset. */
    void add_point_uniform(size_type n) { add_point_uniform(n, dataset_->domain()); }

    void add_point_normal(size_type n, double mean, double stddev) {
        entries_.reserve(n);
        std::normal_distribution<> distr(mean, stddev);
        for (size_type i = 0; i < n; ++i)
            entries_.push_back(new PointQuery<key_type, mapped_type>(std::round(distr(g_))));
    }

    /** Create point queries for the underlying dataset and add them to the workload.
     *  The `n` elements are drawn with a gaussian distribution with the given mean and stddev from the underlying
     *  dataset.
     *  Note: the `mean` and `stddev` correspond to the numbered elements in `dataset_`, i.e. the indices of the elements.
     *  Therefore, the `mean` has to be in range [0, dataset_->size()]. */
    void add_point_normal_index_based(size_type n, double mean, double stddev) {
        entries_.reserve(n);
        std::normal_distribution<> distr(mean, stddev);
        for (size_type i = 0; i < n;) {
            size_type idx = std::round(distr(g_));
            /* Check if idx is inside the data domain. */
            if (not(idx <= dataset_->size())) continue;
            ++i;
            entries_.push_back(new PointQuery<key_type, mapped_type>((*dataset_)[idx].first));
        }
    }

/* ============================================== RANGE ==============================================================*/
    /** Create range queries for the underlying dataset and add them to the workload.
     *  The lower bounds of the `n` elements are drawn uniformly at random in the the domain [dom.min, dom.max -
     *  dataset_size * sel]. The upper bound is chosen according to the given selectivity.
     *  Note: the index_range `range` does *not* refer to the key domain of the dataset but instead to the numbered
     *  elements in `dataset_`.*/
    void add_range_uniform_index_based(size_type n, double sel, index_range dom) {
        assert(dataset_->is_sorted());
        entries_.reserve(n);
        size_type data_size = dataset_->size();
        size_type range_size = sel * data_size;
        assert(range_size <= data_size);
        size_type max = (dom.max < data_size - range_size) ? dom.max : data_size - range_size;
        assert(dom.min <= max);
        std::uniform_int_distribution<size_type> distr(dom.min, max);
        for (size_type i = 0; i < n; ++i) {
            size_type idx = distr(g_);
            entries_.push_back(new RangeQuery<key_type, mapped_type>((*dataset_)[idx].first, (*dataset_)[idx + range_size -1].first));
        }
    }

    /** Create range queries for the underlying dataset and add them to the workload.
     *  The lower bounds of the `n` elements are drawn uniformly at random in the domain [dom.min, dom.max - dataset_size
     *  * sel]. The upper bound is chosen according to the given selectivity. */
    void add_range_uniform_index_based(size_type n, double sel) { add_range_uniform_index_based(n, sel, {0, dataset_->size()}); }

    /** Create range queries for the underlying dataset and add them to the workload.
     *  The lower bounds of the `n` elements are drawn with a gaussian distribution with the given mean and stddev from
     *  the underlying dataset. The upper bound is chosen according to the given selectivity.
     *  Note: the mean and stddev correspond to the numbered elements in `dataset_`, i.e. the indices of the elements. */
    void add_range_normal_index_based(size_type n, double sel, double mean, double stddev) {
        assert(dataset_->is_sorted());
        entries_.reserve(n);
        size_type data_size = dataset_->size();
        size_type range_size = sel * data_size;
        assert(data_size >= (range_size - 1));
        std::normal_distribution<> distr(mean, stddev);
        for (size_type i = 0; i < n;) {
            size_type idx = std::round(distr(g_));
            /* Check if idx is inside the data domain. */
            if (not (idx <= data_size)) continue;
            /* Check if idx + range_size would be outside the data domain. */
            if ((idx + range_size - 1) >= data_size) continue;
            ++i;
            entries_.push_back(new RangeQuery<key_type, mapped_type>((*dataset_)[idx].first, (*dataset_)[idx + range_size - 1].first));
        }
    }

/* ============================================== INSERT =============================================================*/
    /** Create insert queries for the underlying dataset and add them to the workload.
     *  The `n` elements are uniformly distributed in a dense domain beginning at the minimum of the domain `dom`. */
    size_type add_insert_uniform_dense(size_type n, domain_range dom) {
        assert(dataset_->is_sorted());
        entries_.reserve(n);
        size_type insert_count = 0;
        for (auto i = dom.min; i < dom.min + n; i++) {
            /* Create insertion and add to workload. */
            if (std::find(dataset_->begin(), dataset_->end(), std::make_pair(i, i)) != dataset_->end()) {
                /* `key` already contained. */
                continue;
            } else {
                /* Create insertion and add to workload. */
                entries_.push_back(new Insertion<key_type, mapped_type>(i, i));
                ++insert_count;
            }
        }
        has_inserts_ = true;
        return insert_count;
    }

    /** Create unique insert queries for the underlying dataset and add them to the workload.
     *  The `n` elements are drawn uniformly at random in the domain `dom`.
     *  Only keys, which are not already present in the dataset, are added to the workload.
     *  Note: currently, the value is set to the same value as key. */
    size_type add_insert_uniform(size_type n, domain_range dom, bool force=false) {
        assert(dataset_->is_sorted());
        entries_.reserve(n);
        std::uniform_int_distribution<key_type> distr(dom.min, dom.max);
        std::unordered_set<key_type> unique_keys;
        size_type insert_count = 0;
        if (not force) {
            for (size_type i = 0; i < n; ++i) {
                key_type key = distr(g_);
                if (std::find(dataset_->begin(), dataset_->end(), std::make_pair(key, key)) != dataset_->end()) {
                    /* `key` already contained. */
                    continue;
                } else {
                    /* Create insertion and add to workload. */
                    auto unique = unique_keys.insert(key).second; // true if insertion took place, i.e. unique
                    if (unique) {
                        entries_.push_back(new Insertion<key_type, mapped_type>(key, key));
                        ++insert_count;
                    }
                }
            }
        } else {
            /* Retry insert until `n` elements are inserted. */
            while (insert_count < n) {
                key_type key = distr(g_);
                if (std::find(dataset_->begin(), dataset_->end(), std::make_pair(key, key)) != dataset_->end()) {
                    /* `key` already contained. */
                    continue;
                } else {
                    /* Create insertion and add to workload. */
                    auto unique = unique_keys.insert(key).second; // true if insertion took place, i.e. unique
                    if (unique) {
                        entries_.push_back(new Insertion<key_type, mapped_type>(key, key));
                        ++insert_count;
                    }
                }
            }
        }
        has_inserts_ = true;
        return insert_count;
    }

    /** Create insert queries for the underlying dataset and add them to the workload.
     *  The `n` elements are drawn uniformly at random in the data domain.
     *  Only keys, which are not already present in the dataset, are added to the workload. */
    void add_insert_uniform(size_type n, bool force=false) { add_insert_uniform(n, dataset_->domain(), force); }

    /** Create insert queries for the underlying dataset and add them to the workload.
     *  The `n` elements are drawn with a gaussian distribution with the given mean and stddev.
     *  Note: the mean and stddev correspond to actual values in `dataset_`. */
    void add_insert_normal(size_type n, double mean, double stddev) {
        (void) n;
        (void) mean;
        (void) stddev;
        //TODO: implement
        has_inserts_ = true;
    }


/* ============================================== Utils ==============================================================*/
    /* Iterator. */
    using iterator = typename std::vector<wl_type*>::iterator;
    using const_iterator = typename std::vector<wl_type*>::const_iterator;

    iterator begin() { return entries_.begin(); }
    iterator end()   { return entries_.end(); }
    const_iterator cbegin() const { return entries_.cbegin(); }
    const_iterator cend()   const { return entries_.cend(); }

    void shuffle() { std::shuffle(begin(), end(), g_); }

    size_type size() const { return entries_.size(); }
    uint64_t seed() const { return seed_; }

    bool has_inserts() const { return has_inserts_; }

    wl_type & operator[](size_type idx) { return entries_[idx]; }
    const wl_type & operator[](size_type idx) const { return entries_[idx]; }

    /* Splits the workload into train test workload. */
    std::pair<Workload<key_type, mapped_type>*, Workload<key_type, mapped_type>*>
    split_train_test(double train_percentage, bool shuffle=false)
    {
        if (shuffle) this->shuffle();
        Workload<key_type, mapped_type> * train = new Workload<key_type, mapped_type>(dataset_);
        Workload<key_type, mapped_type> * test  = new Workload<key_type, mapped_type>(dataset_);
        size_type cut_idx = train_percentage * size();
        for (size_type i = 0; i < cut_idx; ++i) {
            train->add(entries_[i]);
        }
        for (size_type i = cut_idx; i < size(); ++i) test->add(entries_[i]);
        entries_.clear();
        return std::make_pair(train, test);
    }

/* =========================================== Input/Output ==========================================================*/
    void to_file(const std::string& filename, const std::string& delimiter=":") {
        std::ofstream file;
        file.open(filename);
        if (file.is_open()) {
            file << "Seed" << delimiter << seed_ << std::endl;
            file << "Size" << delimiter << entries_.size() << std::endl;
            file << "Entries" << delimiter << std::endl;
            for(auto e : entries_) file << e->to_string() << "\n";
        } else {
            std::cerr << "Unable to open file." << std::endl;
        }
        file.close();
    }

    static Workload<key_type, mapped_type> * from_file(dataset_type *data, const std::string& filename, const std::string& delimiter=":") {
        std::ifstream file(filename);
        if (file.fail()) throw std::invalid_argument("Workload file does not exist.");
        std::string line;
        if (file.good()) {
            /* Parse first three lines (seed, size, entries). */
            getline(file, line);
            std::string seed = line.substr(0, line.find(delimiter));
            if (seed != "Seed") throw std::invalid_argument("Reading workload from file: First line must specify the seed");
            line.erase(0, line.find(delimiter) + delimiter.length());
            seed = line.substr(0, line.find(delimiter));
            Workload<key_type, mapped_type> * workload = new Workload<key_type, mapped_type>(data, std::stoul(seed));
            getline(file, line);
            std::string size = line.substr(0, line.find(delimiter));
            if (size != "Size") throw std::invalid_argument("Reading workload from file: Second line must specify the size");
            line.erase(0, line.find(delimiter) + delimiter.length());
            size = line.substr(0, line.find(delimiter));
            size_type n = std::stoul(size);
            getline(file, line);
            std::string entries = line.substr(0, line.find(delimiter));
            if (entries != "Entries") throw std::invalid_argument("Reading workload from file: Third line must specify the entries");
            /* Parse workload entries. */
            for (size_type i = 0; i < n; ++i) {
                if (not getline(file, line))
                    throw std::invalid_argument("Reading workload from file: File does not contain the specified amount of entries");
                std::string wl_kind = line.substr(0, line.find(delimiter));
                line.erase(0, line.find(delimiter) + delimiter.length());
                if (wl_kind == "point") {
                    std::string key = line.substr(0, line.find(delimiter));
                    line.erase(0, line.find(delimiter) + delimiter.length());
                    PointQuery<key_type, mapped_type> *p = new PointQuery<key_type, mapped_type>(std::stoul(key));
                    workload->add(p);
                } else if (wl_kind == "range") {
                    std::string lower = line.substr(0, line.find(delimiter));
                    line.erase(0, line.find(delimiter) + delimiter.length());
                    std::string upper = line.substr(0, line.find(delimiter));
                    RangeQuery<key_type, mapped_type> *r = new RangeQuery<key_type, mapped_type>(std::stoul(lower), std::stoul(upper));
                    workload->add(r);
                } else if (wl_kind == "insert") {
                    std::string key = line.substr(0, line.find(delimiter));
                    line.erase(0, line.find(delimiter) + delimiter.length());
                    std::string value = line.substr(0, line.find(delimiter));
                    Insertion<key_type, mapped_type> *in = new Insertion<key_type, mapped_type>(std::stoul(key), std::stoul(value));
                    workload->add(in);
                } else { throw std::invalid_argument("Reading workload form file: Unsupported workload type " + std::string(wl_kind)); }
            }
            return workload;
        } else { throw std::invalid_argument("Reading workload form file: Invalid workload file"); }
    }
};
