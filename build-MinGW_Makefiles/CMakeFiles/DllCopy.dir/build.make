# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = D:\softwares\cmake-3.20.2-windows-x86_64\bin\cmake.exe

# The command to remove a file.
RM = D:\softwares\cmake-3.20.2-windows-x86_64\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\Documents\GitHub\apic2d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\Documents\GitHub\apic2d\build-MinGW_Makefiles

# Utility rule file for DllCopy.

# Include any custom commands dependencies for this target.
include CMakeFiles/DllCopy.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/DllCopy.dir/progress.make

CMakeFiles/DllCopy:
	cd /d D:\Documents\GitHub\apic2d && D:\softwares\cmake-3.20.2-windows-x86_64\bin\cmake.exe -E copy_directory D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-src/bin/ D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/apic2d/Release/

DllCopy: CMakeFiles/DllCopy
DllCopy: CMakeFiles/DllCopy.dir/build.make
.PHONY : DllCopy

# Rule to build all files generated by this target.
CMakeFiles/DllCopy.dir/build: DllCopy
.PHONY : CMakeFiles/DllCopy.dir/build

CMakeFiles/DllCopy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\DllCopy.dir\cmake_clean.cmake
.PHONY : CMakeFiles/DllCopy.dir/clean

CMakeFiles/DllCopy.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\Documents\GitHub\apic2d D:\Documents\GitHub\apic2d D:\Documents\GitHub\apic2d\build-MinGW_Makefiles D:\Documents\GitHub\apic2d\build-MinGW_Makefiles D:\Documents\GitHub\apic2d\build-MinGW_Makefiles\CMakeFiles\DllCopy.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/DllCopy.dir/depend
