# Copyright 2019-2021 Hauke Mewes, Simeon Ehrig
#
# This file is part of vikunja.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 3.18)

# needs to split up the test in two targets, because the nvcc 10.2
# has problems to compile the Catch2 code, which is required for
# the custom main
set(_TARGET_NAME "test_cxx")
set(_LIB_TARGET_NAME "cxxTestCases")

alpaka_add_library(
  ${_LIB_TARGET_NAME}
  src/cxx.cpp
  )
target_link_libraries(${_LIB_TARGET_NAME}
  PRIVATE
  vikunja::internalvikunja
  Catch2::Catch2
  )
target_include_directories(${_LIB_TARGET_NAME}
  PRIVATE
  include)


add_executable(
  ${_TARGET_NAME}
  src/custom_main.cpp
  )
target_link_libraries(${_TARGET_NAME}
  PRIVATE
  Catch2::Catch2
  ${_LIB_TARGET_NAME}
  )
target_include_directories(${_TARGET_NAME}
  PRIVATE
  include)

add_test(NAME ${_TARGET_NAME} COMMAND ${_TARGET_NAME} --cxx ${ALPAKA_CXX_STANDARD} ${_VIKUNJA_TEST_OPTIONS})
