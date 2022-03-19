/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cstring>
#include <iostream>
#include <string>
#include <tudat/utils/details.h>

namespace tudat {
namespace utils {
namespace data {

// verbose level argument
std::string download_file(const char *remote_url,
                          const char *cache,
                          const int verbosity,
                          bool try_unzip,
                          const char *prefix) {

    // prefix_str variable of size MAX_LEN
    char prefix_str[MAX_PATH];
    char filename[MAX_PATH];
    char ret_path[MAX_PATH];

    // if prefix is not specified, use the default, which is homedir/.tudat
    if (prefix == nullptr) {
        // concatenate with /.tudat
        strcpy(prefix_str, details::get_homedir());
        strcat(prefix_str, "/.tudat");

    } else {
        strcpy(prefix_str, prefix);
    }

    // create destination filename
    strcpy(filename, prefix_str);
    strcat(filename, "/");
    strcat(filename, details::url_to_filename(remote_url).c_str());

    if (verbosity > 0) {
        std::cout << "Downloading " << remote_url << " to " << filename
                  << std::endl;
    }

    // if the file is a zip file, then the return path is to the extracted
    // directory
    if (details::is_zip_file(filename) && try_unzip) {
        // return path without .zip
        strcpy(ret_path, details::remove_extension(filename).c_str());
    } else {
        // return path is the determined filename
        strcpy(ret_path, filename);
    }

    // check if cache string == "true"
    if (cache != nullptr && strcmp(cache, "true") == 0) {
        // check if file exists
        if (details::file_exists(filename)) {
            // file exists
            if (verbosity > 0) {
                std::cout << "File cached, skipping download." << std::endl;
            }
            return ret_path;
        } else {
            // file does not exist
            if (verbosity > 0) {
                std::cout << "File not cached, downloading." << std::endl;
            }
            details::download_file_from_source(remote_url, filename);
        }
    } else if (cache != nullptr && strcmp(cache, "update") == 0) {
        // check if file exists
        if (details::file_exists(filename)) {
            // file exists
            if (verbosity > 0) {
                std::cout << "File cached, but forcing update." << std::endl;
            }
        } else {
            // file does not exist
            if (verbosity > 0) {
                std::cout << "File not cached, downloading." << std::endl;
            }
        }
        details::download_file_from_source(remote_url, filename);
    } else if (cache != nullptr && strcmp(cache, "false") == 0) {
        // check if file exists
        if (verbosity > 0) {
            std::cout << "Downloading file" << std::endl;
        }
        details::download_file_from_source(remote_url, filename);
    } else {
        throw std::invalid_argument("Invalid cache argument. Valid arguments are "
                                    "'true', 'update', and 'false'.");
    }

    // if file is zipped, unzip it to the same directory
    if (details::is_zip_file(filename) && try_unzip) {
        details::unzip(filename, ret_path, verbosity);
    }

    // return success
    return ret_path;
}
} // namespace data
} // namespace utils
} // namespace tudat
