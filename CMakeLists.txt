cmake_minimum_required(VERSION 3.10)

project(RastClipper
    VERSION 1.0.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(GNUInstallDirs)

# Library target
add_library(RastClipper
    src/RastClipper.cpp
)

target_include_directories(RastClipper PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Install rules
install(TARGETS RastClipper
    EXPORT RastClipperTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/RastClipper
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT RastClipperTargets
    FILE RastClipperTargets.cmake
    NAMESPACE RastClipper::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RastClipper
)

# Config file for find_package()
include(CMakePackageConfigHelpers)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/RastClipperConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/RastClipperConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RastClipper
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/RastClipperConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/RastClipperConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/RastClipperConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RastClipper
)

# Tests (only when built as top-level project)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    enable_testing()
    add_executable(test_rastclipper tests/test_rastclipper.cpp)
    target_link_libraries(test_rastclipper PRIVATE RastClipper)
    add_test(NAME RastClipperTests COMMAND test_rastclipper)
endif()
