#pragma once

#include <iterator>
#include <map>
#include <memory>
#include <utility>
#include <type_traits>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <typeinfo>

#include "util/hash.hpp"
#include "util/Exceptions.hpp"

#define IMPORT(WHAT) using WHAT = typename Base::WHAT

enum class Node_type { SoA, Tree, Hash };

/** Represents the data structure to store key-value pairs inside a `HorizontalPartition`. */
template<typename Key, typename Value>
struct Node
{
    using size_type = std::size_t;
    using key_type = Key;
    using mapped_type = Value;
    using value_type = std::pair<const key_type, mapped_type>;
    using entry_type = std::pair<key_type, mapped_type>;

    protected:
    Node_type node_type_;

    /* Key Iterator. */
    struct key_iterator_base {
        virtual ~key_iterator_base() { }
        virtual void next() = 0;
        virtual const key_type & get_key() const = 0;
        virtual bool is_equal(key_iterator_base &other) const = 0;
    };

    struct key_iterator {
        private:
        std::unique_ptr<key_iterator_base> impl_;

        public:
        explicit key_iterator (std::unique_ptr<key_iterator_base> &&impl) : impl_(std::move(impl)) { }

        key_iterator & operator++() { impl_->next(); return *this; }
        const key_type & operator*() const { return impl_->get_key(); }
        const key_type & operator->() const { return &impl_->get_key(); }

        bool operator==(key_iterator &other) const { return impl_->is_equal(*other.impl_); }
        bool operator!=(key_iterator &other) const { return not (*this == other); }
    };

    /* Value Iterator. */
    struct value_iterator_base {
        virtual ~value_iterator_base() { }
        virtual void next() = 0;
        virtual mapped_type & get_value() const = 0;
        virtual bool is_equal(value_iterator_base &other) const = 0;
    };

    struct value_iterator {
        private:
        std::unique_ptr<value_iterator_base> impl_;

        public:
        explicit value_iterator (std::unique_ptr<value_iterator_base> &&impl) : impl_(std::move(impl)) { }

        value_iterator & operator++() { impl_->next(); return *this; }
        mapped_type & operator*() const { return impl_->get_value(); }
        mapped_type & operator->() const { return &impl_->get_value(); }

        bool operator==(value_iterator &other) const { return impl_->is_equal(*other.impl_); }
        bool operator!=(value_iterator &other) const { return not (*this == other); }
    };

    public:
    virtual std::optional<mapped_type> get(key_type) = 0;
    virtual bool insert(key_type, mapped_type) = 0;
    virtual bool remove(key_type) = 0;
    virtual bool remove_delete(key_type) = 0;
    virtual size_type size() const = 0;
    virtual size_type capacity() const = 0;
    Node_type node_type() const { return node_type_; }
    bool empty() const { return not size(); }
    virtual ~Node() { };

    virtual key_iterator key_begin() = 0;
    virtual key_iterator key_end() = 0;
    virtual value_iterator value_begin() = 0;
    virtual value_iterator value_end() = 0;

    virtual Node<Key, Value>* clone() = 0;
    virtual Node<Key, Value>* cloneEmpty() = 0;
    virtual std::vector<Key>* get_keys() = 0;
    virtual void clear() = 0;

    virtual void print(std::ostream &out) const = 0;
    virtual void set_capacity(size_type c) = 0;

    friend std::ostream & operator<<(std::ostream &out, Node<Key, Value> &n) {
        n.print(out);
        return out;
    }

    friend std::ostream & operator<<(std::ostream &out, const Node<Key, Value> &n) {
        n.print(out);
        return out;
    }
};

/** Wrapper for the STL container `std::unordered_map`. */
template<typename Key, typename Value>
struct HashNode : public Node<Key, Value>
{
    using Base = Node<Key, Value>;
    using hasher = std::hash<Key>;
    IMPORT(key_iterator_base);
    IMPORT(key_iterator);
    IMPORT(value_iterator_base);
    IMPORT(value_iterator);
    IMPORT(size_type);
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(entry_type);

    private:
    std::unordered_map<key_type, mapped_type, hasher> entries_;
    size_type capacity_ = 4;
    hasher hash_;

    public:
    HashNode(const hasher &hash = hasher()) : hash_(hash) { Base::node_type_ = Node_type::Hash; }
    HashNode(size_type cap, const hasher &hash = hasher())
        : capacity_(cap), hash_(hash)
    {
        Base::node_type_ = Node_type::Hash;
    }

    std::optional<mapped_type> get(key_type key) override {
        auto it = entries_.find(key);
        if (it != entries_.end()) {
            return std::optional<mapped_type>(it->second);
        } else {
            return std::nullopt;
        }
    }

    bool insert(key_type k, mapped_type e) override {
        auto search = entries_.find(k);
        if (search != entries_.end()) {
            /* Element with same key found -> DuplicateKeyException. */
            DuplicateKeyException ex("Insertion of key " + std::to_string(k) + " failed");
            throw ex;
        }
        if (size() >= capacity()) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* Key is not a duplicate and there is space left -> insert. */
        entries_.emplace_hint(search, k, e);
        return true;
    }
    bool remove(key_type k) override { return entries_.erase(k); }

    bool remove_delete(key_type k) override {
        if constexpr (std::is_pointer_v<mapped_type>) {
            // erase and delete
            auto it = entries_.find(k);
            if (it != entries_.end()) {
                delete it->second;
                entries_.erase(it);
                return true;
            }
            return false;
        } else
            return entries_.erase(k);
    }

    std::unordered_map<key_type, mapped_type, hasher> & entries() { return entries_; }

    /* Key Iterator. */
    struct key_iterator_impl : public Node<key_type, mapped_type>::key_iterator_base {
        using iterator_type = typename std::unordered_map<key_type, mapped_type, hasher>::iterator;
        private:
        iterator_type it_;

        public:
        key_iterator_impl(iterator_type it) : it_(it) { }
        void next() override { ++it_; }
        const key_type & get_key() const override { return it_->first; }
        bool is_equal(key_iterator_base &other) const override {
            if (auto key_it = dynamic_cast<key_iterator_impl*>(&other)) return it_ == key_it->it_;
            return false;
        }
    };

    /* Value Iterator. */
    struct value_iterator_impl : public value_iterator_base {
        using iterator_type = typename std::unordered_map<key_type, mapped_type, hasher>::iterator;
        private:
        iterator_type it_;

        public:
        value_iterator_impl(iterator_type it) : it_(it) { }
        void next() override { ++it_; }
        mapped_type & get_value() const override { return it_->second; }
        bool is_equal(value_iterator_base &other) const override {
            if (auto value_it = dynamic_cast<value_iterator_impl*>(&other)) return it_ == value_it->it_;
            return false;
        }
    };

    size_type size() const override { return entries_.size(); }
    size_type capacity() const override { return capacity_; }

    hasher hash_function() const { return hash_; }

    key_iterator key_begin() override { return key_iterator(std::move(std::make_unique<key_iterator_impl>(entries_.begin()))); }
    key_iterator key_end() override { return key_iterator(std::move(std::make_unique<key_iterator_impl>(entries_.end()))); }

    value_iterator value_begin() override { return value_iterator(std::move(std::make_unique<value_iterator_impl>(entries_.begin()))); }
    value_iterator value_end() override { return value_iterator(std::move(std::make_unique<value_iterator_impl>(entries_.end()))); }

    void print(std::ostream &out) const override {
        if constexpr (std::is_pointer_v<mapped_type>) {
            for (auto e : entries_) {
                out << "K:" << e.first << "\n";
                out << *(e.second);
            }
        } else {
            for (auto e : entries_) { out << "K:" << e.first << " - V:" << e.second << "\n"; }
        }
    }

    Node<Key, Value>* cloneEmpty() override {
        return new HashNode<key_type, mapped_type>(this->capacity_);
    }

    std::vector<Key>* get_keys() override {
        std::vector<Key>* keys = new std::vector<Key>();
        for(auto const& elem : entries_) {
            keys->push_back(elem.first);
        }
        return keys;
    }

    Node<key_type, mapped_type>* clone() override {
        HashNode* copy = new HashNode<key_type, mapped_type>(this->capacity_);
        copy->entries_ = std::unordered_map<key_type, mapped_type, hasher>(this->entries_);
        return copy;
    }

    void clear() override { entries_.clear(); }

    void set_capacity(size_type c) override {
        if ( c >= size()) capacity_ = c;
        else throw std::invalid_argument("Capacity reduction beyong current size not allowed");
    }
};

/** Wrapper for the STL container `std::_map` representing a tree structure. */
template<typename Key, typename Value>
struct TreeNode : public Node<Key, Value>
{
    using Base = Node<Key, Value>;
    IMPORT(key_iterator);
    IMPORT(key_iterator_base);
    IMPORT(value_iterator_base);
    IMPORT(value_iterator);
    IMPORT(size_type);
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(entry_type);

    private:
    std::map<key_type, mapped_type> entries_;
    size_type capacity_ = 4;

    public:
    TreeNode() { Base::node_type_ = Node_type::Tree; }
    TreeNode(size_type cap) : capacity_(cap) { Base::node_type_ = Node_type::Tree; }

    std::optional<mapped_type> get(key_type key) override {
        auto it = entries_.find(key);
        if (it != entries_.end()) {
            return std::optional<mapped_type>(it->second);
        } else {
            return std::nullopt;
        }
    }

    bool insert(key_type k, mapped_type e) override {
        auto search = entries_.find(k);
        if (search != entries_.end()) {
            /* Element with same key found -> DuplicateKeyException. */
            DuplicateKeyException ex("Insertion of key " + std::to_string(k) + " failed");
            throw ex;
        }
        if (size() >= capacity()) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* Key is not a duplicate and there is space left -> insert. */
        entries_.insert({k, e});
        return true;
    }

    bool remove(key_type k) override { return entries_.erase(k); }

    bool remove_delete(key_type k) override {
        if constexpr (std::is_pointer_v<mapped_type>) {
            // erase and delete
            auto it = entries_.find(k);
            if (it != entries_.end()) {
                delete it->second;
                entries_.erase(it);
                return true;
            }
            return false;
        } else
            return entries_.erase(k);
    }

    /* Key Iterator. */
    struct key_iterator_impl : public Node<key_type, mapped_type>::key_iterator_base {
        using iterator_type = typename std::map<key_type, mapped_type>::iterator;
        private:
        iterator_type it_;

        public:
        key_iterator_impl(iterator_type it) : it_(it) { }
        void next() override { ++it_; }
        const key_type & get_key() const override { return it_->first; }
        bool is_equal(key_iterator_base &other) const override {
            if (auto key_it = dynamic_cast<key_iterator_impl*>(&other)) return it_ == key_it->it_;
            return false;
        }
    };

    /* Value Iterator. */
    struct value_iterator_impl : public value_iterator_base {
        using iterator_type = typename std::map<key_type, mapped_type>::iterator;
        private:
        iterator_type it_;

        public:
        value_iterator_impl(iterator_type it) : it_(it) { }
        void next() override { ++it_; }
        mapped_type & get_value() const override { return it_->second; }
        bool is_equal(value_iterator_base &other) const override {
            if (auto value_it = dynamic_cast<value_iterator_impl*>(&other)) return it_ == value_it->it_;
            return false;
        }
    };

    std::size_t size() const override { return entries_.size(); }
    size_type capacity() const override { return capacity_; }

    std::map<key_type, mapped_type> & entries() { return entries_; }

    key_iterator key_begin() override { return key_iterator(std::move(std::make_unique<key_iterator_impl>(entries_.begin()))); }
    key_iterator key_end() override { return key_iterator(std::move(std::make_unique<key_iterator_impl>(entries_.end()))); }

    value_iterator value_begin() override { return value_iterator(std::move(std::make_unique<value_iterator_impl>(entries_.begin()))); }
    value_iterator value_end() override { return value_iterator(std::move(std::make_unique<value_iterator_impl>(entries_.end()))); }

    void print(std::ostream &out) const override {
        if constexpr (std::is_pointer_v<mapped_type>) {
            for (auto e : entries_) {
                out << "K:" << e.first << "\n";
                out << *(e.second);
            }
        } else {
            for (auto e : entries_) { out << "K:" << e.first << " - V:" << e.second << "\n"; }
        }
    }

    Node<Key, Value>* cloneEmpty() override { return new TreeNode<key_type, mapped_type>(this->capacity_); }

    std::vector<Key>* get_keys() override {
        std::vector<Key>* keys = new std::vector<Key>();
        for(auto const& elem : entries_) {
            keys->push_back(elem.first);
        }
        return keys;
    }

    Node<key_type, mapped_type> * clone() override {
        //std::map<key_type, mapped_type> entries;
        //for (auto const& [key, value] : this->entries_) { entries.insert({key, value}); }
        TreeNode* copy = new TreeNode(this->capacity_);
        //copy->entries_ = entries;
        copy->entries_ = std::map<key_type, mapped_type>(this->entries_);
        return copy;
    }

    void clear() override { entries_.clear(); }

    void set_capacity(size_type c) override {
        if ( c >= size()) capacity_ = c;
        else throw std::invalid_argument("Capacity reduction beyong current size not allowed");
    }

    std::map<key_type, mapped_type> * get_entries() { return &entries_; }
};

/** Structure of Arrays. */
template<typename Key, typename Value>
struct SoA : public Node<Key, Value>
{
    using Base = Node<Key, Value>;
    IMPORT(key_iterator);
    IMPORT(key_iterator_base);
    IMPORT(value_iterator_base);
    IMPORT(value_iterator);
    IMPORT(size_type);
    IMPORT(key_type);
    IMPORT(mapped_type);
    IMPORT(entry_type);

    private:
    std::vector<key_type> keys_;
    std::vector<mapped_type> values_;
    size_type capacity_ = 4;

    public:
    SoA() { Base::node_type_ = Node_type::SoA; }
    SoA(size_type cap) : capacity_(cap) { Base::node_type_ = Node_type::SoA; }

    std::optional<mapped_type> get(key_type key) override {
        auto k_begin = keys_.begin();
        auto k_end = keys_.end();
        auto k_it = std::find(k_begin, k_end, key);
        if (k_it != k_end) {
            /* Key found. */
            auto distance = std::distance(k_begin, k_it);
            auto e_it = values_.begin();
            std::advance(e_it, distance);
            return std::optional<mapped_type>(*e_it);
        }
        return std::nullopt;
    }

    bool insert(key_type k, mapped_type e) override {
        auto k_begin = keys_.begin();
        auto k_end = keys_.end();
        auto k_it = std::find(k_begin, k_end, k);
        if (k_it != k_end) {
            /* Element with same key found -> DuplicateKeyException. */
            DuplicateKeyException ex("Insertion of key " + std::to_string(k) + " failed");
            throw ex;
        }
        if (size() >= capacity()) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* Key is not a duplicate and there is space left -> insert. */
        keys_.push_back(k);
        values_.push_back(e);
        return true;
    }

    /** Sorted insert. */
    bool insert_sorted(key_type k, mapped_type e) {
        auto k_begin = keys_.begin();
        auto k_end = keys_.end();
        auto k_it = std::lower_bound(k_begin, k_end, k);
        if (k_it != k_end and *k_it == k) {
            DuplicateKeyException ex("Insertion of key " + std::to_string(k) + " failed");
            throw ex;
        }
        if (size() >= capacity()) {
            /* No capacity left -> insertion failed. */
            return false;
        }
        /* Key is not a duplicate and there is space left -> sorted insert. */
        auto k_pos_index = std::distance(k_begin, k_it);
        keys_.emplace(k_it, k);
        auto e_pos = values_.begin();
        std::advance(e_pos, k_pos_index);
        values_.emplace(e_pos, e);
        return true;
    }

    bool remove(key_type k) override {
        auto k_begin = keys_.begin();
        auto k_end = keys_.end();
        auto k_it = std::find(k_begin, k_end, k);
        if (k_it != k_end) {
            auto distance = std::distance(k_begin, k_it);
            auto e_it = values_.begin();
            std::advance(e_it, distance);
            keys_.erase(k_it);
            values_.erase(e_it);
            return true;
        }
        return false;
    }

    bool remove_delete(key_type k) override {
        auto k_begin = keys_.begin();
        auto k_end = keys_.end();
        auto k_it = std::find(k_begin, k_end, k);
        if (k_it != k_end) {
            auto distance = std::distance(k_begin, k_it);
            auto e_it = values_.begin();
            std::advance(e_it, distance);
            keys_.erase(k_it);
            if constexpr (std::is_pointer_v<mapped_type>) {
                delete *e_it;
                values_.erase(e_it);
            } else {
                values_.erase(e_it);
            }
            return true;
        }
        return false;
    }

    /* Key Iterator. */
    struct key_iterator_impl : public key_iterator_base {
        private:
        size_type idx_;
        SoA<key_type, mapped_type> &s_;

        public:
        key_iterator_impl(size_type idx, SoA<key_type, mapped_type> &s) : idx_(idx), s_(s) { }
        void next() override { idx_++; }
        const key_type & get_key() const override { return s_.keys_[idx_]; }
        bool is_equal(key_iterator_base &other) const override {
            if (auto key_it = dynamic_cast<key_iterator_impl*>(&other)) {
                return idx_ == key_it->idx_;
            }
            return false;
        }
    };

    /* Value Iterator. */
    struct value_iterator_impl : public value_iterator_base {
        private:
        size_type idx_;
        SoA<key_type, mapped_type> &s_;

        public:
        value_iterator_impl(size_type idx, SoA<key_type, mapped_type> &s) : idx_(idx), s_(s) { }
        void next() override { idx_++; }
        mapped_type & get_value() const override { return s_.values_[idx_]; }
        bool is_equal(value_iterator_base &other) const override {
            if (auto value_it = dynamic_cast<value_iterator_impl*>(&other)) {
                return idx_ == value_it->idx_;
            }
            return false;
        }
    };

    std::vector<key_type> & keys() { return keys_; }
    std::vector<mapped_type> & values() { return values_; }

    std::size_t size() const override { return keys_.size(); }
    size_type capacity() const override { return capacity_; }

    key_iterator key_begin() override { return key_iterator(std::make_unique<key_iterator_impl>(0, *this)); }
    key_iterator key_end() override { return key_iterator(std::make_unique<key_iterator_impl>(size(), *this)); }

    value_iterator value_begin() override { return value_iterator(std::make_unique<value_iterator_impl>(0, *this)); }
    value_iterator value_end() override { return value_iterator(std::make_unique<value_iterator_impl>(size(), *this)); }

    Node<key_type, mapped_type> * clone() override {
        //std::vector<key_type> keys;
        //for (auto k : this->keys_) { keys.push_back(k); }
        //std::vector<mapped_type> values;
        //for (auto v : this->values_) { values.push_back(v); }
        SoA* copy = new SoA(this->capacity_);
        //copy->keys_ = keys;
        //copy->values_ = values;
        copy->keys_ = std::vector<key_type>(this->keys_);
        copy->values_ = std::vector<mapped_type>(this->values_);
        return copy;
    }

    Node<Key, Value>* cloneEmpty() override { return new SoA<key_type, mapped_type>(this->capacity_); }

    std::vector<Key>* get_keys() override {
        std::vector<Key>* keys = new std::vector<Key>();
        for(auto k : keys_)
            keys->push_back(k);
        return keys;
    }

    void print(std::ostream &out) const override {
        auto keys_it = keys_.begin();
        auto values_it = values_.begin();
        auto keys_end = keys_.end();
        out << "SoA" << std::endl;
        if constexpr (std::is_pointer_v<mapped_type>) {
            for (; keys_it != keys_end; keys_it++, values_it++) {
                out << "K:" << *keys_it << "\n";
                out << *(*values_it);
            }
        } else {
            for (; keys_it != keys_end; keys_it++, values_it++) {
                out << "K:" << *keys_it << " - V:" << *values_it << "\n";
            }
        }
    }

    void clear() override { keys_.clear(); values_.clear(); }

    void set_capacity(size_type c) override { 
        if ( c >= size()) capacity_ = c;
        else throw std::invalid_argument("Capacity reduction beyong current size not allowed");
    }
};
