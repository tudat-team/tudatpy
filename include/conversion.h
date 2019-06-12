//
// Created by ggarrett on 1-6-19.
//

#ifndef TUDATPY_CONVERSION_H
#define TUDATPY_CONVERSION_H

#include <boost/python/dict.hpp>
#include <boost/shared_ptr.hpp>

// [4]
#include <iostream>
#include <vector>
#include <memory>

#include "boost/shared_ptr.hpp"
#include "boost/python.hpp"
#include "boost/python/stl_iterator.hpp"
// [-4]

#include <unordered_map>


// Converts a C++ map to a python dict
// [1] Adapted from https://gist.github.com/octavifs/5362297
template<class K, class V>
boost::python::dict toPythonDict(std::map<K, V> map) {
    typename std::map<K, V>::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}

template<class K, class V>
boost::python::dict toPythonDict(std::unordered_map<K, V> map) {
    typename std::unordered_map<K, V>::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}


// Handles conversions between boost::shared_ptr and std::shared_ptr.
// [2] Sourced from Fozi: https://stackoverflow.com/questions/6326757/conversion-from-boostshared-ptr-to-stdshared-ptr
namespace {
    template<class SharedPointer>
    struct Holder {
        SharedPointer p;

        Holder(const SharedPointer &p) : p(p) {}

        Holder(const Holder &other) : p(other.p) {}

        Holder(Holder &&other) : p(std::move(other.p)) {}

        void operator()(...) { p.reset(); }
    };
}

template<class T>
std::shared_ptr<T> to_std_ptr(const boost::shared_ptr<T> &p) {
    typedef Holder<std::shared_ptr<T>> H;
    if (H *h = boost::get_deleter<H>(p)) {
        return h->p;
    } else {
        return std::shared_ptr<T>(p.get(), Holder<boost::shared_ptr<T>>(p));
    }
}

template<class T>
boost::shared_ptr<T> to_boost_ptr(const std::shared_ptr<T> &p) {
    typedef Holder<boost::shared_ptr<T>> H;
    if (H *h = std::get_deleter<H>(p)) {
        return h->p;
    } else {
        return boost::shared_ptr<T>(p.get(), Holder<std::shared_ptr<T>>(p));
    }
}

// https://stackoverflow.com/questions/13986581/using-boost-python-stdshared-ptr
/* [3] make boost::python understand std::shared_ptr */
namespace boost {
    template<typename T>
    T *get_pointer(std::shared_ptr<T> p) {
        return p.get();
    }
}

// [4] https://gist.github.com/ggarrett13/5733db92c432633302e79bab04dd84ef
namespace tudatpy::conversion {
    template<typename T>
    inline
    std::vector<T> py_list_to_std_vector(const boost::python::object &iterable) {
        return std::vector<T>(boost::python::stl_input_iterator<T>(iterable),
                              boost::python::stl_input_iterator<T>());
    }

    template<class T>
    inline
    boost::python::list std_vector_to_py_list(std::vector<T> vector) {
        typename std::vector<T>::iterator iter;
        boost::python::list list;
        for (iter = vector.begin(); iter != vector.end(); ++iter) {
            list.append(*iter);
        }
        return list;
    }

    template<typename Fn, Fn fn, typename... Args>
    typename std::result_of<Fn(Args...)>::type
    py_list_to_std_vector_wrapper(Args... args) {
        return fn(std::forward<Args>(args)...);
    }

}

// To-python converter for std::map<K,V>
template<typename K, typename V>
struct map_to_python_dict {
    static PyObject *convert(const std::map<K, V> &cppMap) {
        return boost::python::incref(toPythonDict(cppMap).ptr());
    }
};

// To-python converter for std::unordered_map<K,V>
template<typename K, typename V>
struct unordered_map_to_python_dict {
    static PyObject *convert(const std::unordered_map<K, V> &cppUnorderedMap) {
        return boost::python::incref(toPythonDict(cppUnorderedMap).ptr());
    }
};


#endif //TUDATPY_CONVERSION_H
