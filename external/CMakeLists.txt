cmake_minimum_required(VERSION 3.16)

project(Test)

file(DOWNLOAD https://curl.haxx.se/download/curl-7.70.0.tar.gz
        ${CMAKE_CURRENT_SOURCE_DIR}/curl-7.70.0.tar.gz
        STATUS Status
        SHOW_PROGRESS
        )

