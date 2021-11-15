#pragma once

#include <exception>

struct CapacityException : public std::exception
{
    const std::string message;

    CapacityException() : message("Insertion failed due to no capacity left.") {}
    CapacityException(const std::string& m) : message(m) {}

    virtual const char * what() const throw() { return message.c_str(); }
};

struct DuplicateKeyException : public std::exception
{
    const std::string message;

    DuplicateKeyException() : message("Insertion failed due to duplicate key.") {}
    DuplicateKeyException(const std::string& m) : message(m) {}

    virtual const char * what() const throw() { return message.c_str(); }
};

struct SearchMethodNotAvailable : public std::exception
{
    virtual const char * what() const throw() {
        return "IndexSearch strategy not available for given node type.";
    }
};

struct KeyNotFoundException : public std::exception
{
    virtual const char * what() const throw() {
        return "Lookup failed due to missing key.";
    }
};

struct InvalidInsertionException : public std::exception
{
    virtual const char * what() const throw() {
        return "Inserting nested partition not possible";
    }
};

struct SearchMethodNotTrained : public std::exception
{
    virtual const char * what() const throw() {
        return "IndexSearch strategy not trained.";
    }
};

struct FileNotFoundException : public std::exception
{
    virtual const char * what() const throw() {
        return "File not found.";
    }
};
