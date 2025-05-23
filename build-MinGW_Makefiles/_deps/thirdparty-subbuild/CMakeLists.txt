# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.20.2)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(thirdparty-populate NONE)

# Pass through things we've already detected in the main project to avoid
# paying the cost of redetecting them again in ExternalProject_Add()
set(GIT_EXECUTABLE [==[D:/Program Files/Git/mingw64/bin/git.exe]==])
set(GIT_VERSION_STRING [==[2.40.0.windows.1]==])
set_property(GLOBAL PROPERTY _CMAKE_FindGit_GIT_EXECUTABLE_VERSION
  [==[D:/Program Files/Git/mingw64/bin/git.exe;2.40.0.windows.1]==]
)


include(ExternalProject)
ExternalProject_Add(thirdparty-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/nepluno/thirdparty_win64_common.git" "GIT_TAG" "apic2d"
                    SOURCE_DIR          "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-src"
                    BINARY_DIR          "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)


