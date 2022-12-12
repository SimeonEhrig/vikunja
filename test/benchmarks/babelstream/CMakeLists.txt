#
# Copyright 2022 Erik Zenker, Benjamin Worpitz, Jan Stephan, Bernhard Manfred Gruber
#
# This file exemplifies usage of alpaka.
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED “AS IS” AND ISC DISCLAIMS ALL WARRANTIES WITH
# REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL ISC BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
# IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

cmake_minimum_required(VERSION 3.18)

set(_TARGET_NAME_VIKUNJA_BABELSTREAM_ALPAKA "babelstream_alpaka")

alpaka_add_executable(${_TARGET_NAME_VIKUNJA_BABELSTREAM_ALPAKA}
    src/main.cpp
    src/AlpakaStream.cpp)
target_include_directories(${_TARGET_NAME_VIKUNJA_BABELSTREAM_ALPAKA} PRIVATE include)

target_link_libraries(${_TARGET_NAME_VIKUNJA_BABELSTREAM_ALPAKA}
    PRIVATE
    vikunja::internalvikunja
)

target_compile_definitions(${_TARGET_NAME_VIKUNJA_BABELSTREAM_ALPAKA} PUBLIC ALPAKA)

set(_TARGET_NAME_VIKUNJA_BABELSTREAM_VIKUNJA "babelstream_vikunja")

alpaka_add_executable(${_TARGET_NAME_VIKUNJA_BABELSTREAM_VIKUNJA}
    src/main.cpp
    src/VikunjaStream.cpp)
target_include_directories(${_TARGET_NAME_VIKUNJA_BABELSTREAM_VIKUNJA} PRIVATE include)

target_link_libraries(${_TARGET_NAME_VIKUNJA_BABELSTREAM_VIKUNJA}
    PRIVATE
    vikunja::internalvikunja
)

target_compile_definitions(${_TARGET_NAME_VIKUNJA_BABELSTREAM_VIKUNJA} PUBLIC VIKUNJA)
