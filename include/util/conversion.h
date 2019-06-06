//
// Created by ggarrett on 1-6-19.
//

#ifndef TUDATPY_CONVERSION_H
#define TUDATPY_CONVERSION_H

#include <boost/python/dict.hpp>

// Converts a C++ map to a python dict
// Adapted from https://gist.github.com/octavifs/5362297
template <class K, class V>
boost::python::dict toPythonDict(std::map<K, V> map) {
    typename std::map<K, V>::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}

#endif //TUDATPY_CONVERSION_H
