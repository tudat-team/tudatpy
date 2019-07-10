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
#include "typedefs.h"


// Converts a C++ map to a python dict
// [1] Adapted from https://gist.github.com/octavifs/5362297
template<class K, class V>
boost::python::dict toPythonDict(std::map <K, V> map) {
    typename std::map<K, V>::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}

template<class K, class V>
boost::python::dict toPythonDict(std::unordered_map <K, V> map) {
    typename std::unordered_map<K, V>::iterator iter;
    boost::python::dict dictionary;
    for (iter = map.begin(); iter != map.end(); ++iter) {
        dictionary[iter->first] = iter->second;
    }
    return dictionary;
}
//        std::map <std::string, stdBodySettingsPointer> my_map;
//
//        // Both are Python List objects
//        PyObject *pKeys = PyDict_Keys(obj_ptr);
//        PyObject *pValues = PyDict_Values(obj_ptr);
//
//        // Cycles through dictionary.
//        for (Py_ssize_t i = 0; i < PyDict_Size(obj_ptr); ++i) {
//            my_map.insert(std::pair<std::string, stdBodySettingsPointer>(
//                    boost::python::extract<std::string>(PyList_GetItem(pKeys, i)),
//                    boost::python::extract<boost::python::object>(*PyList_GetItem(pValues, i));
//        }
//template<class K, class V>
//std::map fromPythonDict(PyObject * obj_ptr) {
//
//    std::map<K, V> std_map;
//
//    boost::python::dict dictionary;
//
//    for (iter = boost::python::dict.begin(); iter != boost::python::dict.end(); ++iter) {
//        dictionary[iter->first] = iter->second;
//    }

    //        // Both are Python List objects
//        PyObject *pKeys = PyDict_Keys(obj_ptr);
//        PyObject *pValues = PyDict_Values(obj_ptr);
//
//        // Cycles through dictionary.
//        for (Py_ssize_t i = 0; i < PyDict_Size(obj_ptr); ++i) {
//            my_map.insert(std::pair<std::string, stdBodySettingsPointer>(
//                    boost::python::extract<std::string>(PyList_GetItem(pKeys, i)),
//                    boost::python::extract<boost::python::object>(*PyList_GetItem(pValues, i));

//    return dictionary;
//}



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
std::shared_ptr <T> to_std_ptr(const boost::shared_ptr <T> &p) {
    typedef Holder<std::shared_ptr < T>>
    H;
    if (H * h = boost::get_deleter<H>(p)) {
        return h->p;
    } else {
        return std::shared_ptr<T>(p.get(), Holder<boost::shared_ptr < T>>
        (p));
    }
}

template<class T>
boost::shared_ptr <T> to_boost_ptr(const std::shared_ptr <T> &p) {
    typedef Holder<boost::shared_ptr < T>>
    H;
    if (H * h = std::get_deleter<H>(p)) {
        return h->p;
    } else {
        return boost::shared_ptr<T>(p.get(), Holder<std::shared_ptr < T>>
        (p));
    }
}

// https://stackoverflow.com/questions/13986581/using-boost-python-stdshared-ptr
/* [3] make boost::python understand std::shared_ptr */
namespace boost {
    template<typename T>
    T *get_pointer(std::shared_ptr <T> p) {
        return p.get();
    }
}

// [4] https://gist.github.com/ggarrett13/5733db92c432633302e79bab04dd84ef
namespace tudatpy::conversion {
    template<typename T>
    inline
    std::vector <T> py_list_to_std_vector(const boost::python::object &iterable) {
        return std::vector<T>(boost::python::stl_input_iterator<T>(iterable),
                              boost::python::stl_input_iterator<T>());
    }

    template<class T>
    inline
    boost::python::list std_vector_to_py_list(std::vector <T> vector) {
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
    static PyObject *convert(const std::map <K, V> &cppMap) {
        return boost::python::incref(toPythonDict(cppMap).ptr());
    }
};

//std::string py_string_to_std_string(PyASCIIObject *py_file_name)
//{
//    int len = PyUnicode_GET_LENGTH(*py_file_name);     // Not sure how you write that in python.
//    std::string str;
//
//    for(int i = 0; i < len; i++)
//        str[i] = py_file_name[i];
//    return str;
//}

//struct python_bodysettings_dict_to_map {
//    python_bodysettings_dict_to_map() {
//        boost::python::converter::registry::push_back(
//                &convertible,
//                &construct,
//                boost::python::type_id<stdBodySettingsPointerMap>());
//    }
//
//    // Determine if obj_ptr can be converted in to a PyDict
//    static void *convertible(PyObject *obj_ptr) {
//        if (!PyDict_Check(obj_ptr)) return 0;
//        return obj_ptr;
//    }
//
//    // Convert obj_ptr into a QString
//    static void construct(
//            PyObject *obj_ptr,
//            boost::python::converter::rvalue_from_python_stage1_data *data) {
//
//        std::map <std::string, stdBodySettingsPointer> my_map;
//
//        // Both are Python List objects
//        PyObject *pKeys = PyDict_Keys(obj_ptr);
//        PyObject *pValues = PyDict_Values(obj_ptr);
//
//        // Cycles through dictionary.
//        for (Py_ssize_t i = 0; i < PyDict_Size(obj_ptr); ++i) {
//            my_map.insert(std::pair<std::string, stdBodySettingsPointer>(
//                    boost::python::extract<std::string>(PyList_GetItem(pKeys, i)),
//                    boost::python::extract<boost::python::object>(*PyList_GetItem(pValues, i));
//        }
//
//        // Grab pointer to memory into which to construct the new QString
//        void *storage = (
//                (boost::python::converter::rvalue_from_python_storage<stdBodySettingsPointerMap> *)
//                        data)->storage.bytes;
//
//        // in-place construct the new QString using the character data
//        // extraced from the python object
////        new(storage) my_map;
//
//        // Stash the memory chunk pointer for later use by boost.python
//        data->convertible = storage;
//    }
//};
//
// To-python converter for std::unordered_map<K,V>
template<typename K, typename V>
struct unordered_map_to_python_dict {
    static PyObject *convert(const std::unordered_map <K, V> &cppUnorderedMap) {
        return boost::python::incref(toPythonDict(cppUnorderedMap).ptr());
    }
};


//template<typename V>
//struct python_dict_to_unordered_map
//{
//    python_dict_to_unordered_map()
//    {
//        boost::python::converter::registry::push_back(
//                &convertible,
//                &construct,
//                boost::python::type_id());
//    }
//
//    // Determine if obj_ptr can be converted in to a PyDict
//    static void* convertible(PyObject* obj_ptr)
//    {
//        if (!PyDict_Check(obj_ptr)) return 0;
//        return obj_ptr;
//    }
//
//    // Convert obj_ptr into a QString
//    static void construct(
//            PyObject* obj_ptr,
//            boost::python::converter::rvalue_from_python_stage1_data* data)
//    {
//        // Return map c++
//        std::unordered_map<std::string, V> my_map;
//
//        // Python Dictionary object
//        PyObject *pDict = obj_ptr;
//
//        // Both are Python List objects
//        PyObject *pKeys = PyDict_Keys(pDict);
//        PyObject *pValues = PyDict_Values(pDict);
//
//        // Verify that pDict is a dict (should be ensured by convertible())
//        assert(pDict);
//
//        // Cycles through dictionary.
//        for (Py_ssize_t i = 0; i < PyDict_Size(pDict); ++i) {
//
//            // PyString_AsString returns a char*
//            my_map.insert( std::pair<std::string, V>(
//                    *PyString_AsString( PyList_GetItem(pKeys,   i) ),
//                    *PyList_GetItem( pValues, i) );
//        }
//
//        // Grab pointer to memory into which to construct the new QString
//        void* storage = (
//                (boost::python::converter::rvalue_from_python_storage*)
//                        data)->storage.bytes;
//
//        // in-place construct the new QString using the character data
//        // extraced from the python object
//        new (storage) my_map;
//
//        // Stash the memory chunk pointer for later use by boost.python
//        data->convertible = storage;
//    }
//};

// http://cci.lbl.gov/cctbx_sources/scitbx/stl/map_wrapper.h
template <typename MapType>
struct from_python_dict
{
    typedef typename MapType::key_type k_t;
    typedef typename MapType::mapped_type m_t;

    from_python_dict()
    {
        boost::python::converter::registry::push_back(
                &convertible,
                &construct,
                boost::python::type_id<MapType>()
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                , &boost::python::converter::wrap_pytype<&PyDict_Type>::get_pytype
#endif
        );
    }

    static void* convertible(PyObject* obj_ptr)
    {
        return PyDict_Check(obj_ptr) ? obj_ptr : 0;
    }

    static void construct(
            PyObject* obj_ptr,
            boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        boost::python::handle<> obj_hdl(boost::python::borrowed(obj_ptr));
        boost::python::object obj_obj(obj_hdl);
        boost::python::extract<boost::python::dict> obj_proxy(obj_obj);
        boost::python::dict other = obj_proxy();
        void* storage = (
                (boost::python::converter::rvalue_from_python_storage<MapType>*)
                        data)->storage.bytes;
        new (storage) MapType();
        data->convertible = storage;
        MapType& self = *((MapType*)storage);
        boost::python::list keys = other.keys();
        int len_keys = boost::python::len(keys);
        for(int i=0;i<len_keys;i++) {
            boost::python::object key_obj = keys[i];
            boost::python::extract<k_t> key_proxy(key_obj);
            if (!key_proxy.check()) {
                PyErr_SetString(PyExc_KeyError, "Unsuitable type.");
                boost::python::throw_error_already_set();
            }
            boost::python::object value_obj = other[key_obj];
            boost::python::extract<m_t> value_proxy(value_obj);
            if (!value_proxy.check()) {
                PyErr_SetString(PyExc_ValueError, "Unsuitable type.");
                boost::python::throw_error_already_set();
            }
            k_t key = key_proxy();
            m_t value = value_proxy();
            self[key] = value;
        }
    }
};



//struct custom_string_to_python_str
//{
//    static PyObject* convert(custom_string const &s)
//    {
//        return boost::python::incref(boost::python::object(s.value()).ptr());
//    }
//};




#endif //TUDATPY_CONVERSION_H
