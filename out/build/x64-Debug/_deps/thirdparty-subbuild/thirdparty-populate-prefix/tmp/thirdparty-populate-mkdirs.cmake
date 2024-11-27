# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-src"
  "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-build"
  "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-subbuild/thirdparty-populate-prefix"
  "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-subbuild/thirdparty-populate-prefix/tmp"
  "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp"
  "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src"
  "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "D:/Documents/GitHub/apic2d/out/build/x64-Debug/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
