/*    Copyright (c) 2010-2022, Delft University of Technology
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
    snprintf(homedir, MAX_PATH, "%s%s", getenv("HOMEDRIVE"),
             getenv("HOMEPATH"));
#else
    snprintf(homedir, MAX_PATH, "%s", getenv("HOME"));
#endif
    return strdup(homedir);
}

} // namespace details
} // namespace utils
} // namespace tudat

#endif // TUDAT_DETAILS_H
