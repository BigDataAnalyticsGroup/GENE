#pragma once

#include <cassert>
#include <string>

namespace util {

template<typename Key>
struct Range
{
    using key_type = Key;

    key_type min;
    key_type max;

    Range(key_type min, key_type max) : min(min), max(max) { assert(min <= max); }
};


/** https://stackoverflow.com/a/37327537/11734577 */
inline int nth_occurrence(const std::string& str, const std::string& findMe, int nth)
{
    size_t  pos = 0;
    int     cnt = 0;

    while( cnt != nth )
    {
        pos+=1;
        pos = str.find(findMe, pos);
        if ( pos == std::string::npos )
            return -1;
        cnt++;
    }
    return pos;
}

}
