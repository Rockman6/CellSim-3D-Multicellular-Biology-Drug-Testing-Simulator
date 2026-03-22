cmake_minimum_required(VERSION 3.24)
project(CellSim VERSION 1.0 LANGUAGES CXX OBJCXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_OBJCXX_STANDARD 20)

# ══════════════════════════════════════════════════════════════════════════
#  CellSim — Computational Cell Biology Simulator
#  Build configuration for macOS with Metal
# ══════════════════════════════════════════════════════════════════════════

# ── vcpkg integration ────────────────────────────────────────────────────
if(DEFINED ENV{VCPKG_ROOT})
    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
elseif(EXISTS "$ENV{HOME}/vcpkg/scripts/buildsystems/vcpkg.cmake")
    set(CMAKE_TOOLCHAIN_FILE "$ENV{HOME}/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
endif()

# ── Find packages ────────────────────────────────────────────────────────
find_package(glfw3 CONFIG REQUIRED)
find_package(glm CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(implot CONFIG REQUIRED)

# ── macOS Frameworks ─────────────────────────────────────────────────────
find_library(METAL_FRAMEWORK Metal REQUIRED)
find_library(METALKIT_FRAMEWORK MetalKit REQUIRED)
find_library(COCOA_FRAMEWORK Cocoa REQUIRED)
find_library(QUARTZCORE_FRAMEWORK QuartzCore REQUIRED)
find_library(IOKIT_FRAMEWORK IOKit REQUIRED)
find_library(UTI_FRAMEWORK UniformTypeIdentifiers REQUIRED)

# cgltf (header-only GLB parser)
find_path(CGLTF_INCLUDE_DIRS "cgltf.h")

# ── Source files ─────────────────────────────────────────────────────────
set(SOURCES
    src/main.mm
    src/gpu/MetalContext.mm
    src/render/Camera.cpp
    src/render/MeshLibrary.mm
    src/render/GLBLoader.mm
)

# ── Executable ───────────────────────────────────────────────────────────
add_executable(CellSim ${SOURCES})

target_include_directories(CellSim PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CGLTF_INCLUDE_DIRS}
)

target_link_libraries(CellSim PRIVATE
    glfw
    glm::glm
    imgui::imgui
    implot::implot
    ${METAL_FRAMEWORK}
    ${METALKIT_FRAMEWORK}
    ${COCOA_FRAMEWORK}
    ${QUARTZCORE_FRAMEWORK}
    ${IOKIT_FRAMEWORK}
    ${UTI_FRAMEWORK}
)

# ── Compile Metal shaders ────────────────────────────────────────────────
set(METAL_SHADER_DIR "${CMAKE_SOURCE_DIR}/assets/shaders")
set(METAL_OUTPUT_DIR "${CMAKE_BINARY_DIR}")

# Find all .metal files
file(GLOB METAL_SOURCES "${METAL_SHADER_DIR}/*.metal")

# Compile each .metal → .air, then link into default.metallib
set(METAL_AIR_FILES "")
foreach(METAL_SRC ${METAL_SOURCES})
    get_filename_component(METAL_NAME ${METAL_SRC} NAME_WE)
    set(AIR_FILE "${METAL_OUTPUT_DIR}/${METAL_NAME}.air")
    add_custom_command(
        OUTPUT ${AIR_FILE}
        COMMAND xcrun -sdk macosx metal -c ${METAL_SRC} -o ${AIR_FILE}
        DEPENDS ${METAL_SRC}
        COMMENT "Compiling Metal shader: ${METAL_NAME}.metal"
    )
    list(APPEND METAL_AIR_FILES ${AIR_FILE})
endforeach()

add_custom_command(
    OUTPUT "${METAL_OUTPUT_DIR}/default.metallib"
    COMMAND xcrun -sdk macosx metallib ${METAL_AIR_FILES} -o "${METAL_OUTPUT_DIR}/default.metallib"
    DEPENDS ${METAL_AIR_FILES}
    COMMENT "Linking Metal shader library"
)

add_custom_target(MetalShaders ALL
    DEPENDS "${METAL_OUTPUT_DIR}/default.metallib"
)

add_dependencies(CellSim MetalShaders)

# ── Compiler flags ───────────────────────────────────────────────────────
target_compile_options(CellSim PRIVATE
    -Wall -Wextra -Wno-unused-parameter
    $<$<CONFIG:Release>:-O2>
    $<$<CONFIG:Debug>:-O0 -g>
)

# Mark .mm files as Objective-C++
set_source_files_properties(
    src/main.mm
    src/gpu/MetalContext.mm
    src/render/MeshLibrary.mm
    src/render/GLBLoader.mm
    PROPERTIES COMPILE_FLAGS "-x objective-c++"
)

message(STATUS "CellSim — Metal compute + render | Target: ${CMAKE_SYSTEM_PROCESSOR}")
