/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <tudat/utils/details.h>

namespace tudat {
namespace utils {
namespace data {

// verbose level argument
std::string download_file(const char *remote_url, const int verbosity,
                          bool try_unzip, const char *prefix) {

    // prefix_str variable of size MAX_LEN
    char prefix_str[MAX_PATH];
    char filename[MAX_PATH];
    char ret_path[MAX_PATH];
    char curl_cmd[1200];
    char unzip_cmd[1200];
    int uzip_ret;
    int curl_ret;

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

    // create curl command
    strcpy(curl_cmd, "curl -o "); // output to filename
    strcat(curl_cmd, filename);
    strcat(curl_cmd, " -L ");           // allow follow redirect
    strcat(curl_cmd, remote_url);       // remote url
    strcat(curl_cmd, " --create-dirs"); // create directory if not exists
    curl_ret = system(curl_cmd);

    // post download checks
    if (0 != curl_ret) {
        // throw download error
        throw std::runtime_error("Error downloading file from: " +
                                 std::string(remote_url));
    }

    // define ret_path
    if (details::is_zip_file(filename)) {

    } else {
    }

    // if file is zipped, unzip it to the same directory
    if (details::is_zip_file(filename) && try_unzip) {
        // return path without .zip
        strcpy(ret_path, details::remove_extension(filename).c_str());

        strcpy(unzip_cmd, "unzip ");
        if (verbosity > 1) {
            strcat(unzip_cmd, "-v ");
        } else {
            strcat(unzip_cmd, "-q ");
        }
        strcat(unzip_cmd, "-o ");
        strcat(unzip_cmd, filename);
        strcat(unzip_cmd, " -d ");
        strcat(unzip_cmd, ret_path);
        uzip_ret = system(unzip_cmd);
        if (0 != uzip_ret) {
            throw std::runtime_error("Error unzipping file: " +
                                     std::string(filename));
        }
    }

    else {
        // return path without .gz
        strcpy(ret_path, filename);
    }

    // assert that file exists
    if (!details::file_exists(filename)) {
        throw std::runtime_error("Downloaded file does not exist: " +
                                 std::string(filename));
    } else {
        // return success
        return ret_path;
    }
}
} // namespace data
} // namespace utils
} // namespace tudat
