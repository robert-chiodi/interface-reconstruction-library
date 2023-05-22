# Inclusion of GoogleTest
set(IRL_TEST_SOURCE_DIR "${CMAKE_SOURCE_DIR}/tests/src/")


# If path to GoogleTest not given, then download from git
if(NOT GOOGLETEST_DIR)
  # Download and unpack googletest at configure time
  configure_file("${CMAKE_SOURCE_DIR}/cmake/Utils/GoogleTestDownload.cmake" googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()
  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                   ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                   EXCLUDE_FROM_ALL)
else()
  add_subdirectory("${GOOGLETEST_DIR}"
                   "${CMAKE_CURRENT_BINARY_DIR}/googletest-build"
		   EXCLUDE_FROM_ALL)		 		 
endif()


# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
  include_directories(SYSTEM "${gtest_SOURCE_DIR}/include")
endif()

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(irl_test)
target_link_libraries(irl_test irl gtest gmock gtest_main)
target_include_directories(irl_test PRIVATE "${PROJECT_SOURCE_DIR}")
target_include_directories(irl_test SYSTEM PRIVATE "${EIGEN_INCLUDE_DIR}")
include(GoogleTest)
gtest_discover_tests(irl_test)
set_target_properties(irl_test
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/bin"
    COMPILE_FLAGS "${IRL_CXX_FLAGS}"
    CXX_STANDARD 17
)

if(IRL_USE_ABSL)
  target_link_libraries(irl_test absl::inlined_vector absl::flat_hash_map)
  target_include_directories(irl_test SYSTEM PRIVATE "${ABSEIL_DIR}")
endif()


# Add test files to executable. (This is irl_root/tests/src)
add_subdirectory(./src)

