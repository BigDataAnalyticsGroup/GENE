#pragma once

#include <algorithm>
#include <cassert>
#include <cstring>
#include <emmintrin.h> // x86 SSE intrinsics
#include <iostream>
#include <optional>
#include <ostream>
#include <tsl/robin_map.h>
#include <unordered_map>
#include <vector>
#include "util/hash.hpp"

#define IMPORT(WHAT) using WHAT = typename Base::WHAT

namespace opt {

enum class PARTITION_type {
    /* Inner Nodes. */
    InnerSortedBinary,
    InnerSortedLinear,

    /* Leaf Nodes. */
    LeafSortedBinary,
    LeafSortedInterpolation,

    LeafHashUnorderedMap,
    LeafRobinHoodHashTable
};

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct HorizontalPartition
{
    using key_type = Key;
    using mapped_type = T;
    using value_type = std::pair<const key_type, T>;
    using key_compare = Compare;
    using size_type = std::size_t;
    using hp_type = opt::HorizontalPartition<key_type, mapped_type, key_compare>;

    PARTITION_type partition_type;

    protected:
    size_type slotmax_ = 64;
    bool leaf_;

    public:
    virtual ~HorizontalPartition() { }

    size_type capacity() const { return slotmax_; }
    virtual key_type max_key() const = 0;
    bool is_leaf() const { return leaf_; }

    virtual void dot(std::ostream &out, const std::string &id, bool struct_only) const = 0;
};


/*======================================================================================================================
 * Inner Partitions
 *====================================================================================================================*/
template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct InnerPartition : public HorizontalPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);

    InnerPartition() { Base::leaf_ = false; }
    virtual ~InnerPartition() { }

    virtual std::optional<hp_type*> find(const key_type &key) const = 0;
    virtual bool insert(const key_type &key, hp_type *hp) = 0;
};

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct InnerSortedBinary final : public InnerPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);

    private:
    size_type slotuse_ = 0;
    std::vector<key_type> keys_;
    std::vector<hp_type*> partitions_;

    public:
    InnerSortedBinary() = default;
    InnerSortedBinary(size_type slotmax) {
        Base::slotmax_ = slotmax;
        Base::partition_type = PARTITION_type::InnerSortedBinary;
        keys_.reserve(Base::slotmax_);
        partitions_.reserve(Base::slotmax_);
    }

    ~InnerSortedBinary() { for (auto p : partitions_) delete p; }

    private:
    size_type find_lower(const key_type &key) const {
        /* Binary search. */
        if (slotuse_ == 0) return 0;
        size_type lo = 0, hi = slotuse_;
        while (lo < hi) {
            size_type mid = (lo + hi) >> 1;
            if (key <= keys_[mid]) {
                hi = mid; // key <= mid
            } else {
                lo = mid + 1; // key > mid
            }
        }
        return lo;
    }

    public:
    bool insert(const key_type &key, hp_type *hp) {
        // TODO: only use slotmax_ - 1 many keys, fan-out of 64
        if (slotuse_ >= Base::slotmax_) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* There is space left -> sorted insert. */
        auto k_begin = keys_.begin();
        auto k_end = keys_.end();
        auto k_it = std::lower_bound(k_begin, k_end, key);
        auto k_idx = std::distance(k_begin, k_it);
        /* Check if partition would be inserted into the last available slot.  If so , only insert the hp, not the key. 
         * TODO: check if this actually works.
         * Adapt the find method accordingly.
         * Maybe, the right approach is to just implement inserts on leaf nodes and construct inner nodes on demand.
         */
        keys_.emplace(k_it, key);
        auto hp_it = partitions_.begin();
        std::advance(hp_it, k_idx);
        partitions_.emplace(hp_it, hp);
        ++slotuse_;
        return true;
    }

    std::optional<hp_type*> find(const key_type &key) const {
        size_type slot = find_lower(key);
        if (slot == slotuse_) return std::nullopt;
        return std::optional<hp_type*>(partitions_[slot]);
    }

    key_type max_key() const { return keys_[slotuse_ - 1]; } // TODO: fix for empty partition

    void dot(std::ostream &out, const std::string &id, bool struct_only) const {
        /* Write node identifier and setup tabular format. */
        out << "node" << id << "[label=<\n"
            << "<table cellspacing=\"0\" border=\"0\" cellborder=\"1\">\n"
            << "\t<tr>";

        if (struct_only) {
            out << "<td PORT=\"P0\">SortedBSearch</td>";
        } else {
            /* Write keys and label them with the corresponding child partition. */
            for (size_type i = 0; i < keys_.size(); ++i) {
                out << "<td PORT=\"P" << i << "\">" << keys_[i] << "</td>";
            }
        }

        out << "</tr>\n"
            << "</table>>];\n";

        /* Dot recursively. */
        for (size_type i = 0; i < partitions_.size(); ++i) {
            std::string rec_id(id + std::to_string(i));
            partitions_[i]->dot(out, rec_id, struct_only);
            if (struct_only) {
                out << "node" << id << ":P0 -> node" << rec_id << ";\n";
            } else {
                out << "node" << id << ":P" << std::to_string(i) << " -> node" << rec_id << ";\n";
            }
        }
    }
};

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct InnerSortedLinear final : public InnerPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);
    using value_type = std::pair<const key_type, hp_type*>;

    private:
    size_type slotuse_ = 0;
    std::vector<key_type> keys_;
    std::vector<hp_type*> partitions_;

    public:
    InnerSortedLinear() = default;
    InnerSortedLinear(size_type slotmax) {
        Base::slotmax_ = slotmax;
        Base::partition_type = PARTITION_type::InnerSortedLinear;
        keys_.reserve(Base::slotmax_);
        partitions_.reserve(Base::slotmax_);
    }

    ~InnerSortedLinear() { for (auto p : partitions_) delete p; }

    private:
    size_type find_lower(const key_type &key) const {
        /* Linear search. */
        for (auto i = 0; i < slotuse_; ++i) if (key <= keys_[i]) return i;
        return slotuse_;
    }

    public:
    bool insert(const key_type &key, hp_type *hp) {
        if (slotuse_ >= Base::slotmax_) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* There is space left -> sorted insert. */
        auto k_begin = keys_.begin();
        auto k_idx = find_lower(key); // linear search
        keys_.emplace(k_begin + k_idx, key);
        auto v_it = partitions_.begin();
        std::advance(v_it, k_idx);
        partitions_.emplace(v_it, hp);
        ++slotuse_;
        return true;
    }

    std::optional<hp_type*> find(const key_type &key) const {
        size_type slot = find_lower(key);
        if (slot == slotuse_) return std::nullopt;
        return std::optional<hp_type*>(partitions_[slot]);
    }

    key_type max_key() const { return keys_[slotuse_ - 1]; } // TODO: fix for empty partition

    value_type operator[](size_type idx) { return {keys_[idx], partitions_[idx]}; }

    void dot(std::ostream &out, const std::string &id, bool struct_only) const {
        /* Write node identifier and setup tabular format. */
        out << "node" << id << "[label=<\n"
            << "<table cellspacing=\"0\" border=\"0\" cellborder=\"1\">\n"
            << "\t<tr>";

        if (struct_only) {
            out << "<td PORT=\"P0\">SortedLinearSearch</td>";
        } else {
            /* Write keys and label them with the corresponding child partition. */
            for (size_type i = 0; i < keys_.size(); ++i) {
                out << "<td PORT=\"P" << i << "\">" << keys_[i] << "</td>";
            }
        }

        out << "</tr>\n"
            << "</table>>];\n";

        /* Dot recursively. */
        for (size_type i = 0; i < partitions_.size(); ++i) {
            std::string rec_id(id + std::to_string(i));
            partitions_[i]->dot(out, rec_id, struct_only);
            if (struct_only) {
                out << "node" << id << ":P0 -> node" << rec_id << ";\n";
            } else {
                out << "node" << id << ":P" << std::to_string(i) << " -> node" << rec_id << ";\n";
            }
        }
    }
};


/*======================================================================================================================
 * Leaf Partitions
 *====================================================================================================================*/
template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct LeafPartition : public HorizontalPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);

    private:
    LeafPartition<key_type, mapped_type, key_compare> *next_ = nullptr;

    public:
    LeafPartition() { Base::leaf_ = true; }
    virtual ~LeafPartition() { }

    virtual std::optional<mapped_type> find(const key_type &key) const = 0;
    /** Returns a std::pair consisting of:
     *  first:  the number of keys in the range [lower, upper]
     *  second: true iff we do *not* have to scan additional leaves after this one
     */
    virtual std::pair<size_type, bool> find(const key_type &lower, const key_type &upper) const = 0;
    virtual bool insert(const key_type &key, const mapped_type &value) = 0;

    LeafPartition<key_type, mapped_type, key_compare> * get_next() const { return next_; }
    void set_next(LeafPartition<key_type, mapped_type, key_compare> *next) { next_ = next; }
};

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct LeafSortedBinary final : public LeafPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);

    private:
    size_type slotuse_ = 0;
    std::vector<key_type> keys_;
    std::vector<mapped_type> values_;

    public:
    LeafSortedBinary() = default;
    LeafSortedBinary(size_type slotmax) {
        Base::slotmax_ = slotmax;
        Base::partition_type = PARTITION_type::LeafSortedBinary;
        keys_.reserve(Base::slotmax_);
        values_.reserve(Base::slotmax_);
    }

    private:
    size_type find_lower(const key_type &key) const {
        /* Binary search. */
        if (slotuse_ == 0) return 0;
        size_type lo = 0, hi = slotuse_;
        while (lo < hi) {
            size_type mid = (lo + hi) >> 1;
            if (key <= keys_[mid]) {
                hi = mid; // key <= mid
            } else {
                lo = mid + 1; // key > mid
            }
        }
        return lo;
    }

    public:
    bool insert(const key_type &key, const mapped_type &v) {
        if (slotuse_ >= Base::slotmax_) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* There is space left -> sorted insert. */
        auto k_begin = keys_.begin();
        auto k_end = keys_.end();
        auto k_it = std::lower_bound(k_begin, k_end, key);
        keys_.emplace(k_it, key);
        auto k_idx = std::distance(k_begin, k_it);
        auto v_it = values_.begin();
        std::advance(v_it, k_idx);
        values_.emplace(v_it, v);
        ++slotuse_;
        return true;
    }

    std::optional<mapped_type> find(const key_type &key) const {
        size_type slot = find_lower(key);
        return (slot < slotuse_ && keys_[slot] == key)
            ? std::optional<mapped_type>(values_[slot]) : std::nullopt;
    }

    std::optional<mapped_type> lower_bound(const key_type &key) const {
        size_type slot = find_lower(key);
        return (slot < slotuse_) ? std::optional<mapped_type>(values_[slot]) : std::nullopt;
    }

    std::pair<size_type, bool> find(const key_type &lower, const key_type &upper) const {
        size_type res = 0;
        size_type i = find_lower(lower); // set `i` to first element not less than `lower`
        for (auto end = slotuse_; i < end; ++i) {
            /* We know that keys_[i] >= lower due to while loop. */
            if (keys_[i] <= upper) ++res;
            else return {res, true}; // early abort, we do not need to scan the next leaf
        }
        return {res, false};
    }

    key_type max_key() const { return keys_[slotuse_ - 1]; } // TODO: fix for empty partition

    void dot(std::ostream &out, const std::string &id, bool struct_only) const {
        /* Write node identifier and setup tabular format. */
        out << "node" << id << "[label=<\n"
            << "<table cellspacing=\"0\" border=\"0\" cellborder=\"1\">\n";

        if (struct_only) {
            out << "\t<tr><td>SortedBSearch</td>";
        } else {
            /* Write keys. */
            out << "\t<tr>";
            for (size_type i = 0; i < keys_.size(); ++i) {
                out << "<td>" << keys_[i] << "</td>";
            }
            out << "</tr>\n";

            /* Write values. */
            out << "\t<tr>";
            for (size_type i = 0; i < values_.size(); ++i) {
                out << "<td>" << values_[i] << "</td>";
            }
        }

        out << "</tr>\n"
            << "</table>>];\n";
    }
};

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct LeafSortedInterpolation final : public LeafPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(value_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);

    private:
    size_type slotuse_ = 0;
    std::vector<key_type> keys_;
    std::vector<mapped_type> values_;

    public:
    LeafSortedInterpolation() = default;
    LeafSortedInterpolation(size_type slotmax) {
        Base::slotmax_ = slotmax;
        Base::partition_type = PARTITION_type::LeafSortedInterpolation;
        keys_.reserve(Base::slotmax_);
        values_.reserve(Base::slotmax_);
    }

    private:
    size_type find_lower(const key_type &key) const {
        /* Interpolation search. */
        auto slotuse_ = keys_.size();
        if (slotuse_ == 0) return 0;
        size_type lo = 0, hi = slotuse_ - 1;
        key_type key_left = keys_[lo];
        key_type key_right = keys_[hi];

        while (lo < hi and key >= key_left and key <= key_right) {
            double slope = (double) (hi - lo) / (key_right - key_left);
            size_type expected = lo + (key - key_left) * slope;
            if (key <= keys_[expected]) {
                hi = expected - 1;
                key_right = keys_[hi];
            } else {
                lo = expected + 1;
                key_left = keys_[lo];
            }
        }
        if (key <= key_left) return lo;
        return hi + 1;
    }

    public:
    bool insert(const key_type &key, const mapped_type &v) {
        if (slotuse_ >= Base::slotmax_) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* There is space left -> sorted insert. */
        auto k_begin = keys_.begin();
        auto k_idx = find_lower(key); // interpolation search
        keys_.emplace(k_begin + k_idx, key);
        auto v_it = values_.begin();
        std::advance(v_it, k_idx);
        values_.emplace(v_it, v);
        ++slotuse_;
        return true;
    }

    std::optional<mapped_type> find(const key_type &key) const {
        size_type slot = find_lower(key);
        return (slot < slotuse_ && keys_[slot] == key)
            ? std::optional<mapped_type>(values_[slot]) : std::nullopt;
    }

    std::optional<mapped_type> lower_bound(const key_type &key) const {
        size_type slot = find_lower(key);
        return (slot < slotuse_) ? std::optional<mapped_type>(values_[slot]) : std::nullopt;
    }

    std::pair<size_type, bool> find(const key_type &lower, const key_type &upper) const {
        /* TODO: Implement. */
        return {0, false};
    }

    value_type operator[](size_type idx) { return {keys_[idx], values_[idx]}; }
    //const value_type & operator[](size_type idx) const { return {keys_[idx], values_[idx]}; } // TODO: does not work

    key_type max_key() const { return keys_[slotuse_ - 1]; } // TODO: fix for empty partition

    void dot(std::ostream &out, const std::string &id, bool struct_only) const {
        /* Write node identifier and setup tabular format. */
        out << "node" << id << "[label=<\n"
            << "<table cellspacing=\"0\" border=\"0\" cellborder=\"1\">\n";

        if (struct_only) {
            out << "\t<tr><td>SortedIntrplSearch</td>";
        } else {
            /* Write keys. */
            out << "\t<tr>";
            for (size_type i = 0; i < keys_.size(); ++i) {
                out << "<td>" << keys_[i] << "</td>";
            }
            out << "</tr>\n";

            /* Write values. */
            out << "\t<tr>";
            for (size_type i = 0; i < values_.size(); ++i) {
                out << "<td>" << values_[i] << "</td>";
            }
        }
        out << "</tr>\n"
            << "</table>>];\n";
    }
};

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct LeafHashUnorderedMap final : public LeafPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    using hasher = Murmur3;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);

    private:
    hasher hasher_;
    size_type slotuse_ = 0;
    key_type max_key_ = 0;
    std::unordered_map<key_type, mapped_type, hasher> entries_;

    public:
    LeafHashUnorderedMap() = default;
    LeafHashUnorderedMap(size_type slotmax) {
        Base::slotmax_ = slotmax;
        Base::partition_type = PARTITION_type::LeafHashUnorderedMap;
        entries_.reserve(Base::slotmax_);
    }

    std::optional<mapped_type> find(const key_type &key) const {
        auto search = entries_.find(key);
        if (search != entries_.end()) return std::optional<mapped_type>(search->second);
        return std::nullopt;
    }
    std::pair<size_type, bool> find(const key_type &lower, const key_type &upper) const {
        auto res = 0;
        for (auto it = entries_.begin(), end = entries_.end(); it != end; ++it) {
            auto k = it->first;
            if (k >= lower && k <= upper) ++res;
        }
        if (upper < max_key_) return {res, true};
        return {res, false};
    }
    bool insert(const key_type &key, const mapped_type &value) {
        if (slotuse_ >= Base::slotmax_) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        if (key > max_key_) max_key_ = key;
        ++slotuse_;
        return entries_.insert({key, value}).second;
    }

    key_type max_key() const { return max_key_; } // TODO: fix for empty partition

    hasher hash_function() const { return hasher_; }

    void dot(std::ostream &out, const std::string &id, bool struct_only) const {
        /* Write node identifier and setup tabular format. */
        out << "node" << id << "[label=<\n"
            << "<table cellspacing=\"0\" border=\"0\" cellborder=\"1\">\n";

        if (struct_only) {
            out << "\t<tr><td>HTUnorderedMap</td>";
        } else {
            /* Write keys. */
            out << "\t<tr>";
            for (auto it = entries_.begin(), end = entries_.end(); it != end; ++it) {
                out << "<td>" << it->first << "</td>";
            }
            out << "</tr>\n";

            /* Write values. */
            out << "\t<tr>";
            for (auto it = entries_.begin(), end = entries_.end(); it != end; ++it) {
                out << "<td>" << it->second << "</td>";
            }
        }

        out << "</tr>\n"
            << "</table>>];\n";
    }
};

template<typename Key,
         typename T,
         typename Compare = std::less<Key>>
struct LeafRobinHoodHashTable final : public LeafPartition<Key, T, Compare>
{
    using Base = HorizontalPartition<Key, T, Compare>;
    using hasher = Murmur3;
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(key_compare);
    IMPORT(size_type);
    IMPORT(hp_type);

    private:
    hasher hasher_;
    size_type slotuse_ = 0;
    key_type max_key_ = 0;
    tsl::robin_map<key_type, mapped_type, hasher> entries_;

    public:
    LeafRobinHoodHashTable() = default;
    LeafRobinHoodHashTable(size_type slotmax) {
        Base::slotmax_ = slotmax;
        Base::partition_type = PARTITION_type::LeafRobinHoodHashTable;
        entries_.reserve(Base::slotmax_);
    }

    std::optional<mapped_type> find(const key_type &key) const {
        auto search = entries_.find(key);
        if (search != entries_.end()) return std::optional<mapped_type>(search->second);
        return std::nullopt;
    }
    std::pair<size_type, bool> find(const key_type &lower, const key_type &upper) const {
        auto res = 0;
        for (auto it = entries_.begin(), end = entries_.end(); it != end; ++it) {
            auto k = it->first;
            if (k >= lower && k <= upper) ++res;
        }
        if (upper < max_key_) return {res, true};
        return {res, false};
    }
    bool insert(const key_type &key, const mapped_type &value) {
        if (slotuse_ >= Base::slotmax_) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        if (key > max_key_) max_key_ = key;
        ++slotuse_;
        return entries_.insert({key, value}).second;
    }

    key_type max_key() const { return max_key_; } // TODO: fix for empty partition

    hasher hash_function() const { return hasher_; }

    void dot(std::ostream &out, const std::string &id, bool struct_only) const {
        /* Write node identifier and setup tabular format. */
        out << "node" << id << "[label=<\n"
            << "<table cellspacing=\"0\" border=\"0\" cellborder=\"1\">\n";

        if (struct_only) {
            out << "\t<tr><td>HTRobinHood</td>";
        } else {
            /* Write keys. */
            out << "\t<tr>";
            for (auto it = entries_.begin(), end = entries_.end(); it != end; ++it) {
                out << "<td>" << it->first << "</td>";
            }
            out << "</tr>\n";

            /* Write values. */
            out << "\t<tr>";
            for (auto it = entries_.begin(), end = entries_.end(); it != end; ++it) {
                out << "<td>" << it->second << "</td>";
            }
        }

        out << "</tr>\n"
            << "</table>>];\n";
    }
};
#undef IMPORT
}
