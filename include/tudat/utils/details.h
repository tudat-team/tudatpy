/*    Copyright (c) 2010-2020, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DETAILS_H
#define TUDAT_DETAILS_H

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>

#define MAX_PATH 255

namespace tudat {

namespace utils {
namespace details {

std::string remove_extension(const std::string &filename) {
    size_t lastdot = filename.find_last_of(".");
    if (lastdot == std::string::npos)
        return filename;
    return filename.substr(0, lastdot);
}

bool file_exists(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

void download_file_from_source(const char *remote_url, const char *filename) {
    char curl_cmd[1200];
    int curl_ret;

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
    } else if (!file_exists(filename)) {
        throw std::runtime_error("Downloaded file does not exist: " +
                                 std::string(filename));
    }
}

void unzip(const char *zip_file, const char *output_dir, int verbosity) {
    char unzip_cmd[1200];
    int unzip_ret;

    // create unzip command
    strcpy(unzip_cmd, "unzip ");
    if (verbosity > 1) {
        strcat(unzip_cmd, "-v ");
    } else {
        strcat(unzip_cmd, "-q ");
    }
    strcat(unzip_cmd, "-o ");
    strcat(unzip_cmd, output_dir);
    strcat(unzip_cmd, " -d ");
    strcat(unzip_cmd, output_dir);
    unzip_ret = system(unzip_cmd);

    // post unzip checks
    if (0 != unzip_ret) {
        // throw unzip error
        throw std::runtime_error("Error unzipping file: " +
                                 std::string(zip_file));
    } else if (!file_exists(output_dir)) {
        throw std::runtime_error("Unzipped directory does not exist: " +
                                 std::string(output_dir));
    }
}

bool is_zip_file(const char *fileName) {
    std::ifstream infile(fileName, std::ios::binary);
    if (!infile.good()) {
        return false;
    }

    char buffer[2];
    infile.read(buffer, 2);
    infile.close();

    return (buffer[0] == 'P' && buffer[1] == 'K');
}

std::string url_to_filename(const std::string &remote_url) {
    std::string file_name = remote_url;
    std::string::size_type last_slash_pos = file_name.rfind('/');
    if (last_slash_pos != std::string::npos) {
        file_name = file_name.substr(last_slash_pos + 1);
    }

    return file_name;
}

// https://cboard.cprogramming.com/c-programming/164689-how-get-users-home-directory.html
static inline char *get_homedir(void) {
    char homedir[MAX_PATH];
#ifdef _WIN32
    snprintf(
        homedir, MAX_PATH, "%s%s", getenv("HOMEDRIVE"), getenv("HOMEPATH"));
#else
    snprintf(homedir, MAX_PATH, "%s", getenv("HOME"));
#endif
    return strdup(homedir);
}

} // namespace details
} // namespace utils
} // namespace tudat

#endif // TUDAT_DETAILS_H
