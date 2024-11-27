
if(NOT "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp/thirdparty-populate-gitinfo.txt" IS_NEWER_THAN "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp/thirdparty-populate-gitclone-lastrun.txt")
  message(STATUS "Avoiding repeated git clone, stamp file is up to date: 'D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp/thirdparty-populate-gitclone-lastrun.txt'")
  return()
endif()

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to remove directory: 'D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-src'")
endif()

# try the clone 3 times in case there is an odd git clone issue
set(error_code 1)
set(number_of_tries 0)
while(error_code AND number_of_tries LESS 3)
  execute_process(
    COMMAND "D:/Program Files/Git/mingw64/bin/git.exe"  clone --no-checkout --config "advice.detachedHead=false" "https://github.com/nepluno/thirdparty_win64_common.git" "thirdparty-src"
    WORKING_DIRECTORY "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps"
    RESULT_VARIABLE error_code
    )
  math(EXPR number_of_tries "${number_of_tries} + 1")
endwhile()
if(number_of_tries GREATER 1)
  message(STATUS "Had to git clone more than once:
          ${number_of_tries} times.")
endif()
if(error_code)
  message(FATAL_ERROR "Failed to clone repository: 'https://github.com/nepluno/thirdparty_win64_common.git'")
endif()

execute_process(
  COMMAND "D:/Program Files/Git/mingw64/bin/git.exe"  checkout apic2d --
  WORKING_DIRECTORY "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-src"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to checkout tag: 'apic2d'")
endif()

set(init_submodules TRUE)
if(init_submodules)
  execute_process(
    COMMAND "D:/Program Files/Git/mingw64/bin/git.exe"  submodule update --recursive --init 
    WORKING_DIRECTORY "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-src"
    RESULT_VARIABLE error_code
    )
endif()
if(error_code)
  message(FATAL_ERROR "Failed to update submodules in: 'D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-src'")
endif()

# Complete success, update the script-last-run stamp file:
#
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy
    "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp/thirdparty-populate-gitinfo.txt"
    "D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp/thirdparty-populate-gitclone-lastrun.txt"
  RESULT_VARIABLE error_code
  )
if(error_code)
  message(FATAL_ERROR "Failed to copy script-last-run stamp file: 'D:/Documents/GitHub/apic2d/build-MinGW_Makefiles/_deps/thirdparty-subbuild/thirdparty-populate-prefix/src/thirdparty-populate-stamp/thirdparty-populate-gitclone-lastrun.txt'")
endif()

